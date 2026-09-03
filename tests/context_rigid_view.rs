//! Contract tests for `src/context/rigid_view.rs`
//! (`.claude/specs/stage-pipeline-02-view.md`).
//!
//! `RigidView` is the single home of the rigid degrees of freedom: the flat
//! `6 * nmol` placement vector (COM block first, Euler block second, three
//! values per molecule each) plus the three operations that cross between it
//! and a `PackContext` —
//!
//! - `write_xcart` — `xcart = com + R(euler) · coor` (the outbound rebuild,
//!   inherited verbatim from `initial::init_xcart_from_x`),
//! - `install_seed` — inject a seed placement + conformer as plain data,
//! - `capture_from_xcart` — read lab-frame coordinates back into
//!   `(com, euler = 0)` plus a centered conformer in `ctx.coor` (the growth
//!   writeback contract).
//!
//! Fixtures build a `PackContext` directly (`PackContext::new` + the public
//! layout fields), the same way `tests/geometry_cache.rs` does — no engine,
//! no solver, no growth driver. Categories: basics, edge cases, immutability
//! (a fresh view is all zeros; a clone is an independent snapshot), plus one
//! hard-coded regression scenario. No physics is asserted here: the view is
//! pure bookkeeping over a coordinate layout.

use molpack::context::RigidView as ContextRigidView;
use molpack::euler::{compcart, eulerrmat};
use molpack::{F, PackContext, RigidView};

/// The crate-root re-export and the `context` path name one type, not two.
/// Compile-time only.
fn _both_paths_are_one_type(v: RigidView) -> ContextRigidView {
    v
}

// ── fixtures ───────────────────────────────────────────────────────────────

/// Two three-atom copies of one type: `ntotat = 6`, `ntotmol = 2`, layout
/// `idfirst = [0]`. `coor` holds one reference conformer **per copy**,
/// sharing `xcart`'s index space (the `PackContext` convention).
fn two_copies_of_three() -> PackContext {
    let mut ctx = PackContext::new(6, 2, 1);
    ctx.nmols = vec![2];
    ctx.natoms = vec![3];
    ctx.idfirst = vec![0];
    ctx.coor = vec![[0.0; 3]; 6];
    ctx
}

/// The regression conformer: two copies whose reference blocks are already
/// centred (each block sums to exactly zero), so the COM survives the
/// `write_xcart` → `capture_from_xcart` round trip.
const CENTERED_COOR: [[F; 3]; 6] = [
    [1.0, 0.0, -0.5],
    [-0.5, 0.25, 0.75],
    [-0.5, -0.25, -0.25],
    [0.5, -1.5, 0.25],
    [-1.25, 0.75, -1.0],
    [0.75, 0.75, 0.75],
];

/// Centroid of `coords` computed with the growth writeback's own arithmetic
/// (`src/grow/driver.rs`): accumulate component-wise in atom order, then
/// divide once by the atom count. Same order, same bits.
fn centroid(coords: &[[F; 3]]) -> [F; 3] {
    let mut com = [0.0 as F; 3];
    for p in coords {
        for k in 0..3 {
            com[k] += p[k];
        }
    }
    for v in com.iter_mut() {
        *v /= coords.len() as F;
    }
    com
}

fn assert_close(got: [F; 3], want: [F; 3], tol: F, what: &str) {
    for k in 0..3 {
        assert!(
            (got[k] - want[k]).abs() <= tol,
            "{what}: component {k} is {} but expected {} (tol {tol})",
            got[k],
            want[k]
        );
    }
}

// ── Category: basics — the flat-vector layout contract ─────────────────────

/// `fresh(nmol)` replaces `vec![0.0; 6 * nmol]`: the view owns its buffer,
/// it is `6 * nmol` long, and every degree of freedom starts at zero.
#[test]
fn rigid_view_fresh_is_zeroed_and_sized() {
    let view = RigidView::fresh(4);
    assert_eq!(view.nmol(), 4, "fresh(4) addresses 4 molecules");
    assert_eq!(
        view.as_slice().len(),
        6 * 4,
        "the flat vector holds 6 variables per molecule (3 COM + 3 Euler)"
    );
    assert!(
        view.as_slice().iter().all(|&v| v == 0.0),
        "a fresh view is all zeros: {:?}",
        view.as_slice()
    );
    for i in 0..4 {
        assert_eq!(view.com(i), [0.0; 3], "fresh com({i})");
        assert_eq!(view.euler(i), [0.0; 3], "fresh euler({i})");
    }
}

/// The view writes the Packmol flat-vector convention: COM of molecule `i` at
/// `x[3*i .. 3*i+3]`, Euler of molecule `i` at
/// `x[3*nmol + 3*i .. 3*nmol + 3*i + 3]`.
///
/// Ported from `tests/grow.rs::placements_view_layout`, which pinned the same
/// contract on the deleted `PlacementsMut`.
#[test]
fn rigid_view_layout() {
    let nmol = 3;
    let mut view = RigidView::fresh(nmol);
    assert_eq!(view.nmol(), nmol);
    view.set_com(1, [1.0, 2.0, 3.0]);
    view.set_euler(2, [0.1, 0.2, 0.3]);
    // Read back through the typed accessors.
    assert_eq!(view.com(1), [1.0, 2.0, 3.0]);
    assert_eq!(view.euler(2), [0.1, 0.2, 0.3]);
    // Untouched copies read zero.
    assert_eq!(view.com(0), [0.0; 3]);
    assert_eq!(view.euler(0), [0.0; 3]);

    // Raw slots: COM block first, Euler block at 3*nmol.
    let x = view.as_slice();
    assert_eq!(&x[3..6], &[1.0, 2.0, 3.0], "com(1) slot");
    assert_eq!(
        &x[15..18],
        &[0.1, 0.2, 0.3],
        "euler(2) slot = x[3*3 + 3*2 ..]"
    );
    // Every other slot untouched.
    for (i, &v) in x.iter().enumerate() {
        if !(3..6).contains(&i) && !(15..18).contains(&i) {
            assert_eq!(v, 0.0, "slot {i} must be untouched");
        }
    }
}

/// `as_mut_slice` is the raw door the GENCAN phases drive; writes through it
/// are visible to the typed accessors and vice versa — one buffer, not two.
#[test]
fn rigid_view_as_mut_slice_shares_one_buffer() {
    let mut view = RigidView::fresh(2);
    view.as_mut_slice()[4] = 7.5; // com(1).y
    view.as_mut_slice()[6] = -1.25; // euler(0).beta
    assert_eq!(view.com(1), [0.0, 7.5, 0.0]);
    assert_eq!(view.euler(0), [-1.25, 0.0, 0.0]);

    view.set_com(0, [1.0, 2.0, 3.0]);
    assert_eq!(&view.as_slice()[0..3], &[1.0, 2.0, 3.0]);
    assert_eq!(view.as_mut_slice().len(), 12);
}

// ── Category: basics — write_xcart (com + R(euler) · coor) ─────────────────

/// `write_xcart` expands the rigid DOF into lab-frame coordinates with the
/// crate's own Euler convention: for every copy, `xcart = com + R · coor`,
/// enumerated over the `(idfirst, natoms, nmols)` layout.
///
/// The expectation is recomputed here from `euler::eulerrmat` + `compcart` —
/// the same arithmetic in the same order, so equality is exact.
#[test]
fn rigid_view_write_xcart_composes_com_and_rotation() {
    let mut ctx = two_copies_of_three();
    ctx.coor = CENTERED_COOR.to_vec();

    let mut view = RigidView::fresh(2);
    let coms = [[10.0, -3.0, 2.5], [-4.0, 7.25, 0.5]];
    let eulers = [[0.3, -0.7, 1.1], [2.0, 0.5, -1.25]];
    for i in 0..2 {
        view.set_com(i, coms[i]);
        view.set_euler(i, eulers[i]);
    }

    view.write_xcart(&mut ctx);

    for imol in 0..2 {
        let (v1, v2, v3) = eulerrmat(eulers[imol][0], eulers[imol][1], eulers[imol][2]);
        for a in 0..3 {
            let icart = 3 * imol + a;
            let want = compcart(&coms[imol], &CENTERED_COOR[icart], &v1, &v2, &v3);
            for (k, (got, want)) in ctx.xcart[icart].iter().zip(&want).enumerate() {
                assert_eq!(
                    got, want,
                    "atom {icart} component {k}: xcart must be com + R(euler)·coor"
                );
            }
        }
    }
}

/// Zero Euler angles are the identity rotation, so `write_xcart` degenerates
/// to a pure translation of the stored conformer. This is the exact shape the
/// growth path relies on: `capture_from_xcart` leaves `euler = 0`, and the
/// outbound rebuild must then reproduce the captured coordinates.
#[test]
fn rigid_view_write_xcart_identity_rotation_is_translation() {
    let mut ctx = two_copies_of_three();
    ctx.coor = CENTERED_COOR.to_vec();

    let mut view = RigidView::fresh(2);
    view.set_com(0, [1.0, 2.0, 4.0]);
    view.set_com(1, [-8.0, 0.5, 16.0]);

    view.write_xcart(&mut ctx);

    for imol in 0..2 {
        let com = view.com(imol);
        for a in 0..3 {
            let icart = 3 * imol + a;
            let want = [
                com[0] + CENTERED_COOR[icart][0],
                com[1] + CENTERED_COOR[icart][1],
                com[2] + CENTERED_COOR[icart][2],
            ];
            for (k, (got, want)) in ctx.xcart[icart].iter().zip(&want).enumerate() {
                assert_eq!(
                    got, want,
                    "atom {icart} component {k}: euler = 0 must be a pure translation"
                );
            }
        }
    }
}

// ── Category: basics — install_seed (plain-data injection) ─────────────────

/// `install_seed` is the constructing form: the seed arrives as two bare
/// slices (never as an `entry`-layer type), the conformer is copied into the
/// leading `coor` block bitwise, and the returned view carries the seed
/// placements verbatim. Atoms past the seed's length (fixed targets) keep
/// whatever `coor` held.
#[test]
fn rigid_view_install_seed_copies_conformer_and_placements() {
    let mut ctx = PackContext::new(8, 2, 1);
    ctx.nmols = vec![2];
    ctx.natoms = vec![3];
    ctx.idfirst = vec![0];
    // Six free-atom slots plus two trailing slots that must not be touched.
    ctx.coor = vec![[9.0, 9.0, 9.0]; 8];

    let seed_coor: [[F; 3]; 6] = CENTERED_COOR;
    let seed_x: [F; 12] = [
        1.5, -2.5, 3.5, // com of molecule 0
        -0.25, 0.75, 8.0, // com of molecule 1
        0.1, 0.2, 0.3, // euler of molecule 0
        -0.4, 0.5, -0.6, // euler of molecule 1
    ];

    let view = RigidView::install_seed(&seed_x, &seed_coor, &mut ctx);

    assert_eq!(view.nmol(), 2, "the seed covers both free copies");
    assert_eq!(
        view.as_slice(),
        &seed_x[..],
        "the seed placements are injected verbatim (zero-conversion chaining)"
    );
    for (i, want) in seed_coor.iter().enumerate() {
        for (k, (got, want)) in ctx.coor[i].iter().zip(want).enumerate() {
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "coor[{i}] component {k} must be a bitwise copy of the seed"
            );
        }
    }
    for i in 6..8 {
        assert_eq!(
            ctx.coor[i],
            [9.0, 9.0, 9.0],
            "coor[{i}] is past the seed and must be left alone"
        );
    }
    // The typed accessors read the injected layout, not a re-derivation.
    assert_eq!(view.com(1), [-0.25, 0.75, 8.0]);
    assert_eq!(view.euler(0), [0.1, 0.2, 0.3]);
}

// ── Category: basics — capture_from_xcart (the growth writeback) ───────────

/// The writeback contract, lifted from the growth drivers: per copy the COM
/// is the centroid of that copy's lab-frame atoms, `ctx.coor` receives the
/// centered conformer (`xcart − com`), and the Euler angles are zeroed —
/// growth produces no rigid rotation, the whole shape lives in `coor`.
///
/// Coordinates are deliberately non-representable in binary, so the assertion
/// also pins the summation order (accumulate in atom order, divide once).
#[test]
fn rigid_view_capture_from_xcart_centers_each_copy() {
    let mut ctx = two_copies_of_three();
    let xcart: [[F; 3]; 6] = [
        [0.1, 0.2, 0.3],
        [1.3, -0.7, 2.9],
        [-0.4, 3.1, 0.55],
        [7.7, 1.1, -2.2],
        [8.3, 0.9, -1.05],
        [9.15, 2.35, -3.4],
    ];
    ctx.xcart = xcart.to_vec();
    // Pre-fill `coor` with a sentinel so the writeback is visible.
    ctx.coor = vec![[42.0; 3]; 6];

    let mut view = RigidView::fresh(2);
    // A stale rotation must be cleared, not preserved.
    view.set_euler(0, [1.0, 2.0, 3.0]);
    view.set_euler(1, [-1.0, -2.0, -3.0]);

    view.capture_from_xcart(&mut ctx);

    for imol in 0..2 {
        let block = &xcart[3 * imol..3 * imol + 3];
        let com = centroid(block);
        let got_com = view.com(imol);
        for (k, (got, want)) in got_com.iter().zip(&com).enumerate() {
            assert_eq!(
                got, want,
                "molecule {imol} COM component {k} must be the centroid of its atoms"
            );
        }
        assert_eq!(
            view.euler(imol),
            [0.0; 3],
            "molecule {imol}: growth writes no rigid rotation, so euler is zeroed"
        );
        for a in 0..3 {
            let icart = 3 * imol + a;
            for k in 0..3 {
                assert_eq!(
                    ctx.coor[icart][k],
                    xcart[icart][k] - com[k],
                    "coor[{icart}] component {k} must be the centered conformer"
                );
            }
        }
    }
    assert_eq!(
        ctx.xcart,
        xcart.to_vec(),
        "capture reads xcart; it must not rewrite it"
    );
}

// ── Category: edge cases ───────────────────────────────────────────────────

/// `fresh(0)`: an empty view is legal (a run with no free molecules), it has
/// no slots, and both context crossings are no-ops.
#[test]
fn rigid_view_fresh_zero_molecules() {
    let mut view = RigidView::fresh(0);
    assert_eq!(view.nmol(), 0);
    assert!(view.as_slice().is_empty(), "fresh(0) has no slots");

    let mut ctx = PackContext::new(0, 0, 0);
    ctx.nmols = Vec::new();
    ctx.natoms = Vec::new();
    ctx.idfirst = Vec::new();
    ctx.coor = Vec::new();

    view.write_xcart(&mut ctx);
    assert!(ctx.xcart.is_empty(), "nothing to expand");
    view.capture_from_xcart(&mut ctx);
    assert!(ctx.coor.is_empty(), "nothing to capture");
    assert_eq!(view.nmol(), 0);
}

/// A molecule index outside `0..nmol` is a programming error, not a runtime
/// condition, and the view must say so.
///
/// Ported from `tests/grow.rs::placements_view_rejects_bad_len`: the deleted
/// `PlacementsMut` could be handed a mis-sized backing slice, while
/// `RigidView::fresh` owns its buffer (illegal state unrepresentable), so the
/// only reachable length error left is an out-of-range molecule.
///
/// This needs an explicit `i < nmol` check on the accessors. The offset
/// arithmetic alone does NOT catch it: `set_com(nmol, ..)` lands on
/// `x[3*nmol .. 3*nmol+3]`, which is molecule 0's **Euler** slot — a
/// silent cross-block write, the worst possible outcome (`set_euler` out of
/// range does trap, since it runs off the end of the buffer; that asymmetry
/// is exactly why the check has to be explicit).
#[test]
#[should_panic]
fn rigid_view_set_com_out_of_range_panics() {
    let mut view = RigidView::fresh(2);
    view.set_com(2, [1.0, 2.0, 3.0]);
}

/// The mirror of the COM case on the Euler block: an out-of-range molecule is
/// refused, never wrapped around into another molecule's slot.
#[test]
#[should_panic]
fn rigid_view_set_euler_out_of_range_panics() {
    let mut view = RigidView::fresh(2);
    view.set_euler(2, [0.1, 0.2, 0.3]);
}

/// A one-atom copy has its COM exactly on the atom and a conformer of exactly
/// zero — no drift from the centroid division.
#[test]
fn rigid_view_capture_from_xcart_single_atom_copy() {
    let mut ctx = PackContext::new(2, 2, 1);
    ctx.nmols = vec![2];
    ctx.natoms = vec![1];
    ctx.idfirst = vec![0];
    ctx.coor = vec![[7.0; 3]; 2];
    ctx.xcart = vec![[1.25, -3.5, 0.125], [-9.0, 0.0, 4.75]];

    let mut view = RigidView::fresh(2);
    view.capture_from_xcart(&mut ctx);

    assert_eq!(view.com(0), [1.25, -3.5, 0.125]);
    assert_eq!(view.com(1), [-9.0, 0.0, 4.75]);
    assert_eq!(ctx.coor[0], [0.0; 3], "a single atom sits on its own COM");
    assert_eq!(ctx.coor[1], [0.0; 3]);
    assert_eq!(view.euler(0), [0.0; 3]);
    assert_eq!(view.euler(1), [0.0; 3]);
}

/// With several types of different sizes, the view enumerates the context's
/// `(idfirst, natoms, nmols)` layout, so a molecule index in `x` is exactly
/// the molecule's index in the `xcart` layout — type-major, copy-major.
#[test]
fn rigid_view_capture_from_xcart_multitype_layout() {
    // type 0: two copies of 2 atoms (icart 0..4); type 1: one copy of 3
    // atoms (icart 4..7).
    let mut ctx = PackContext::new(7, 3, 2);
    ctx.nmols = vec![2, 1];
    ctx.natoms = vec![2, 3];
    ctx.idfirst = vec![0, 4];
    ctx.coor = vec![[0.0; 3]; 7];
    ctx.xcart = vec![
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0], // molecule 0 → centroid [1, 0, 0]
        [0.0, 10.0, 0.0],
        [0.0, 14.0, 0.0], // molecule 1 → centroid [0, 12, 0]
        [0.0, 0.0, 30.0],
        [3.0, 0.0, 30.0],
        [0.0, 3.0, 30.0], // molecule 2 → centroid [1, 1, 30]
    ];

    let mut view = RigidView::fresh(3);
    view.capture_from_xcart(&mut ctx);

    assert_eq!(view.com(0), [1.0, 0.0, 0.0], "molecule 0 = type 0, copy 0");
    assert_eq!(view.com(1), [0.0, 12.0, 0.0], "molecule 1 = type 0, copy 1");
    assert_eq!(view.com(2), [1.0, 1.0, 30.0], "molecule 2 = type 1, copy 0");
    // The centered conformer lands in each copy's OWN `coor` block.
    assert_eq!(ctx.coor[2], [0.0, -2.0, 0.0], "copy 1's first atom");
    assert_eq!(ctx.coor[3], [0.0, 2.0, 0.0], "copy 1's second atom");
    assert_eq!(ctx.coor[4], [-1.0, -1.0, 0.0], "type 1's first atom");
    for i in 0..3 {
        assert_eq!(view.euler(i), [0.0; 3], "molecule {i} euler");
    }
}

// ── Category: immutability ─────────────────────────────────────────────────

/// The view owns its buffer, so a clone is an independent snapshot: writing
/// through the clone leaves the original untouched. (`Debug` is part of the
/// pinned surface too.)
#[test]
fn rigid_view_clone_is_independent_snapshot() {
    let mut original = RigidView::fresh(2);
    original.set_com(0, [1.0, 2.0, 3.0]);
    original.set_euler(1, [0.4, 0.5, 0.6]);

    let mut clone = original.clone();
    assert_eq!(
        clone.as_slice(),
        original.as_slice(),
        "a clone starts equal to its source"
    );

    clone.set_com(0, [-1.0, -2.0, -3.0]);
    assert_eq!(
        original.com(0),
        [1.0, 2.0, 3.0],
        "writing through the clone must not reach the original"
    );
    assert_eq!(clone.com(0), [-1.0, -2.0, -3.0]);
    assert_eq!(
        original.euler(1),
        [0.4, 0.5, 0.6],
        "the untouched half of the original is unchanged"
    );
    assert!(!format!("{original:?}").is_empty(), "RigidView is Debug");
}

// ── Regression scenario (hard-coded golden) ────────────────────────────────

/// Hard-coded golden for the two context crossings, at 1e-12.
///
/// Provenance: the `XCART_GOLDEN` literals were captured on 2026-09-03 from
/// the build at commit c8fb40e — before `RigidView` existed — by running
/// today's `molpack::initial::init_xcart_from_x` on exactly this fixture in a
/// scratch integration test and printing each component with `{:?}` (Rust's
/// shortest round-trip form). `write_xcart` inherits that function verbatim,
/// so the numbers must not move.
///
/// The reverse leg cannot round-trip the Euler angles — `capture_from_xcart`
/// reports a rotation-free copy whose whole shape has been folded into
/// `coor`, which is the writeback contract, not a loss. What must hold is:
/// the COM comes back (the reference blocks are centred), the Euler angles
/// are exactly zero, and rebuilding from the captured state reproduces the
/// same lab-frame coordinates. The `coor` round trip itself is asserted on
/// the zero-rotation leg below, where `R` is the identity.
#[test]
fn rigid_view_regression_xcart_and_capture_golden() {
    const TOL: F = 1e-12;
    const XCART_GOLDEN: [[F; 3]; 6] = [
        [11.10410275417801, -2.8278964923973584, 2.5365717225106734],
        [9.147598602963882, -2.6147818338858686, 2.4956614718464545],
        [9.748298642858106, -3.5573216737167725, 2.467766805642872],
        [-5.467684077808935, 7.842387310633623, 0.7397513752753073],
        [-2.36539646517135, 7.73462779568072, 0.9671265177436698],
        [
            -4.1669194570197154,
            6.1729848936856575,
            -0.20687789301897705,
        ],
    ];
    let coms = [[10.0, -3.0, 2.5], [-4.0, 7.25, 0.5]];
    let eulers = [[0.3, -0.7, 1.1], [2.0, 0.5, -1.25]];

    // ── forward: (com, euler, coor) → xcart ────────────────────────────────
    let mut ctx = two_copies_of_three();
    ctx.coor = CENTERED_COOR.to_vec();
    let mut view = RigidView::fresh(2);
    for i in 0..2 {
        view.set_com(i, coms[i]);
        view.set_euler(i, eulers[i]);
    }
    view.write_xcart(&mut ctx);
    for (i, want) in XCART_GOLDEN.iter().enumerate() {
        assert_close(ctx.xcart[i], *want, TOL, &format!("xcart[{i}]"));
    }

    // ── reverse: xcart → (com, euler = 0, centered coor) ───────────────────
    view.capture_from_xcart(&mut ctx);
    for (i, want) in coms.iter().enumerate() {
        assert_close(view.com(i), *want, TOL, &format!("recaptured com({i})"));
        assert_eq!(
            view.euler(i),
            [0.0; 3],
            "capture reports a rotation-free copy"
        );
    }

    // ── rebuild from the captured state reproduces the same coordinates ────
    view.write_xcart(&mut ctx);
    for (i, want) in XCART_GOLDEN.iter().enumerate() {
        assert_close(ctx.xcart[i], *want, TOL, &format!("rebuilt xcart[{i}]"));
    }

    // ── zero-rotation leg: here `coor` itself round-trips ──────────────────
    let mut ctx0 = two_copies_of_three();
    ctx0.coor = CENTERED_COOR.to_vec();
    let mut view0 = RigidView::fresh(2);
    for (i, com) in coms.iter().enumerate() {
        view0.set_com(i, *com);
    }
    view0.write_xcart(&mut ctx0);
    view0.capture_from_xcart(&mut ctx0);
    for (i, want) in coms.iter().enumerate() {
        assert_close(view0.com(i), *want, TOL, &format!("identity com({i})"));
    }
    for (i, want) in CENTERED_COOR.iter().enumerate() {
        assert_close(ctx0.coor[i], *want, TOL, &format!("identity coor[{i}]"));
    }
}
