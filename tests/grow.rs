//! Integration tests for the chain-growth solver spec
//! (`.claude/specs/chain-growth-solver.md`), running on the `CbmcGrow`
//! entry (`.claude/specs/engine-entry-split.md`).
//!
//! ── Section: Task 1 — named rejections + result surface ───────────────────
//!
//! Covers the entry seam: the named `GrowError` rejections (no silent
//! degradation, spec principle 3) and `PackResult::softened`. NO growth
//! algorithm is exercised here. The rigid-placement layout contract lives in
//! `tests/context_rigid_view.rs` (`rigid_view_layout` /
//! `rigid_view_set_com_out_of_range_panics`).
//!
//! Later sections: Task 2 (internal-coordinate round-trips + random-vars
//! invariants), Task 3 (overlap field), Task 4 (torsion priors / C∞
//! calibration — RED until `src/grow/prior.rs` grows `sample` +
//! `three_state_from_c_inf`). Growth driver + chain statistics +
//! determinism land with Tasks 5-10.

use molpack::grow::field::{OverlapField, Probe};
use molpack::grow::internal::InternalTree;
use molpack::grow::{GrowConfig, GrowError, TorsionPrior};
use molpack::{
    CbmcGrow, F, GenCanPack, Handler, InsideSphereRestraint, PackContext, PackEngine, PackError,
    PackResult, StepInfo, Target, TopologyError,
};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use ndarray::Array1;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use std::sync::{Arc, Mutex};

/// Planar zigzag bead-chain coordinates: tetrahedral (109.5°) bond angles in
/// the x–z plane, all torsions trans. The zigzag is load-bearing: a collinear
/// chain makes every torsion a no-op.
fn zigzag_coords(n: usize, bond_len: F) -> Vec<[F; 3]> {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    (0..n)
        .map(|i| [i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }])
        .collect()
}

/// The `(i, i+1)` bond list of a linear chain.
fn chain_bonds(n: usize) -> Vec<(u32, u32)> {
    (0..n as u32 - 1).map(|i| (i, i + 1)).collect()
}

/// Coordinates + explicit bond list as a `molrs::Frame`, following the
/// `chain_frame` pattern in `tests/optimizer.rs`: atoms block with x/y/z
/// columns and (unless `bonds` is empty) a bonds block with atomi/atomj.
/// Deliberately NO `bond_type` column — this is the shape a PDB CONECT list
/// or a hand-built coarse-grain frame has.
fn frame_from_parts(coords: &[[F; 3]], bonds: &[(u32, u32)]) -> Frame {
    let mut atoms = Block::new();
    for (name, k) in [("x", 0), ("y", 1), ("z", 2)] {
        let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
        atoms
            .insert(name, Array1::from_vec(col).into_dyn())
            .expect("coordinate column");
    }

    let mut frame = Frame::new();
    frame.insert("atoms", atoms);

    if !bonds.is_empty() {
        let mut block = Block::new();
        let ai: Vec<u32> = bonds.iter().map(|&(i, _)| i).collect();
        let aj: Vec<u32> = bonds.iter().map(|&(_, j)| j).collect();
        block
            .insert("atomi", Array1::from_vec(ai).into_dyn())
            .expect("atomi column");
        block
            .insert("atomj", Array1::from_vec(aj).into_dyn())
            .expect("atomj column");
        frame.insert("bonds", block);
    }

    frame
}

/// Zigzag bead chain as a `molrs::Frame` (see [`zigzag_coords`] /
/// [`frame_from_parts`]).
fn chain_frame(n: usize, bond_len: F, with_bonds: bool) -> Frame {
    let bonds = if with_bonds {
        chain_bonds(n)
    } else {
        Vec::new()
    };
    frame_from_parts(&zigzag_coords(n, bond_len), &bonds)
}

/// A generous periodic box that trivially fits two 5-bead chains.
const BOX_MAX: [F; 3] = [20.0, 20.0, 20.0];

// ── CbmcGrow — named rejections (no silent degradation) ───────────────────

/// `Grow` on a `from_coords`-built target (template == None) is a named
/// error, not a silent fall-back to rigid-body packing.
#[test]
fn grow_rejects_target_without_template() {
    // The frozen prior surface: all three variants construct, and config /
    // prior are Debug + Clone.
    let priors = [
        TorsionPrior::Uniform,
        TorsionPrior::Template { kappa: 2.0 },
        TorsionPrior::States(vec![(0.0, 0.645), (2.094, 0.1775), (-2.094, 0.1775)]),
    ];
    let cfg = GrowConfig::new(priors[0].clone());
    let _ = format!("{:?}", cfg.clone());

    let coords: [[F; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.5, 0.0, 1.1],
        [4.0, 0.0, 1.1],
    ];
    let target = Target::from_coords(&coords, &[1.0; 4], 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow without a template must be rejected");

    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::MissingTemplate,
                ..
            }
        ),
        "expected PackError::Grow(MissingTemplate), got {err:?}"
    );
}

/// `Grow` on a template frame with an atoms block but NO bonds block is a
/// named `NoBonds` error — chemistry comes from the bond graph, and molpack
/// does not guess it.
#[test]
fn grow_rejects_template_without_bonds() {
    let target = Target::new(chain_frame(5, 1.5, false), 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow on a bond-less template must be rejected");

    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::Topology(TopologyError::NoBonds),
                ..
            }
        ),
        "expected PackError::Grow(NoBonds), got {err:?}"
    );
}

/// `Grow` on a bonded template with fewer than 3 atoms is a named
/// `TemplateTooSmall` error carrying the actual atom count.
#[test]
fn grow_rejects_tiny_template() {
    let target = Target::new(chain_frame(2, 1.5, true), 2);

    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(42)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("Grow on a 2-atom template must be rejected");

    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::TemplateTooSmall(2),
                ..
            }
        ),
        "expected PackError::Grow(TemplateTooSmall(2)), got {err:?}"
    );
}

// ── Section: Task 2 — InternalTree: internal-coordinate decomposition ──────
//
// Review tests for `src/grow/internal.rs` (spec Design §4a): the round-trip
// (template torsion values must rebuild the template coordinates exactly),
// the random-vars invariant (any free-variable values must preserve every
// bond length and bonded angle — the only detector for a ring-closure bond
// misclassified as a free variable), and the 1-`exclusion_depth` table.

/// Row-major identity rotation: `place_seed` under it, with `origin` at the
/// template position of `seed_atoms()[0]`, reproduces the template's own seed
/// geometry.
const IDENTITY: [[F; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn vsub(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vdot(a: [F; 3], b: [F; 3]) -> F {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vdist(a: [F; 3], b: [F; 3]) -> F {
    let d = vsub(a, b);
    vdot(d, d).sqrt()
}

/// Bonded angle `a–vertex–c` in radians.
fn vangle(a: [F; 3], vertex: [F; 3], c: [F; 3]) -> F {
    let u = vsub(a, vertex);
    let v = vsub(c, vertex);
    (vdot(u, v) / (vdot(u, u).sqrt() * vdot(v, v).sqrt()))
        .clamp(-1.0, 1.0)
        .acos()
}

fn max_abs_dev(a: &[[F; 3]], b: &[[F; 3]]) -> F {
    a.iter()
        .zip(b)
        .flat_map(|(p, q)| (0..3).map(move |k| (p[k] - q[k]).abs()))
        .fold(0.0, F::max)
}

/// Every bonded angle triple `(i, j, k)` — two bonds sharing vertex `j`.
fn bonded_angle_triples(n: usize, bonds: &[(u32, u32)]) -> Vec<(usize, usize, usize)> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(i, j) in bonds {
        adj[i as usize].push(j as usize);
        adj[j as usize].push(i as usize);
    }
    let mut out = Vec::new();
    for (j, nbrs) in adj.iter().enumerate() {
        for a in 0..nbrs.len() {
            for b in a + 1..nbrs.len() {
                out.push((nbrs[a], j, nbrs[b]));
            }
        }
    }
    out
}

/// Rebuild the full coordinate set: seed at the template's own seed position
/// under the identity rotation, then every step under `vars`.
fn rebuild_coords(tree: &InternalTree, template: &[[F; 3]], vars: &[F]) -> Vec<[F; 3]> {
    let mut coords = vec![[F::NAN; 3]; tree.n_atoms()];
    tree.place_seed(template[tree.seed_atoms()[0]], &IDENTITY, &mut coords);
    for k in 0..tree.n_steps() {
        tree.place_step(k, vars, &mut coords);
    }
    coords
}

/// The template's own value for every free variable, collected step by step
/// through the public accessors (`step_var` / `template_var`).
fn template_vars(tree: &InternalTree) -> Vec<F> {
    let mut vars = vec![0.0 as F; tree.n_vars()];
    for k in 0..tree.n_steps() {
        if let Some(v) = tree.step_var(k) {
            vars[v] = tree.template_var(k);
        }
    }
    vars
}

/// Decompose, rebuild with the template's own variable values, and assert the
/// exact round-trip plus the step-partition consistency (seed + steps place
/// every atom exactly once). Returns the tree for extra assertions.
fn assert_roundtrip(coords: &[[F; 3]], bonds: &[(u32, u32)], what: &str) -> InternalTree {
    let tree = InternalTree::from_frame(&frame_from_parts(coords, bonds))
        .unwrap_or_else(|e| panic!("{what}: template must decompose, got {e}"));
    assert_eq!(tree.n_atoms(), coords.len(), "{what}: n_atoms");

    let mut all: Vec<usize> = tree.seed_atoms().to_vec();
    for k in 0..tree.n_steps() {
        let atoms: Vec<usize> = tree.step_atoms(k).collect();
        assert_eq!(
            atoms.len(),
            tree.step_len(k),
            "{what}: step_len({k}) must match step_atoms({k})"
        );
        all.extend(atoms);
    }
    all.sort_unstable();
    assert_eq!(
        all,
        (0..tree.n_atoms()).collect::<Vec<_>>(),
        "{what}: seed + step atoms must cover every atom exactly once"
    );

    let rebuilt = rebuild_coords(&tree, coords, &template_vars(&tree));
    let dev = max_abs_dev(coords, &rebuilt);
    assert!(
        dev < 1e-9,
        "{what}: template-vars round-trip ‖Δ‖∞ = {dev:e}, must be < 1e-9"
    );
    tree
}

/// 10-atom zigzag backbone with two one-atom side branches at tetrahedral-ish
/// positions off backbone atoms 3 and 6. Branch bonds are terminal, so the
/// non-terminal backbone bonds stay the only free variables.
fn branched_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let mut coords = zigzag_coords(10, 1.53);
    let mut bonds = chain_bonds(10);
    let c3 = coords[3];
    coords.push([c3[0], c3[1] + 1.44, c3[2] + 0.51]); // atom 10
    bonds.push((3, 10));
    let c6 = coords[6];
    coords.push([c6[0], c6[1] - 1.44, c6[2] - 0.51]); // atom 11
    bonds.push((6, 11));
    (coords, bonds)
}

/// Chair-ish (non-planar) 6-ring with ~1.54 Å bonds plus a 4-atom zigzag tail
/// off ring atom 0. Ring bonds are cyclic and must never become free
/// variables; the three non-terminal tail bonds must.
fn ring_tail_parts() -> (Vec<[F; 3]>, Vec<(u32, u32)>) {
    let pi = std::f64::consts::PI as F;
    let (r, h) = (1.46 as F, 0.25 as F);
    let mut coords: Vec<[F; 3]> = (0..6)
        .map(|i| {
            let phi = i as F * pi / 3.0;
            [
                r * phi.cos(),
                r * phi.sin(),
                if i % 2 == 0 { h } else { -h },
            ]
        })
        .collect();
    let mut bonds: Vec<(u32, u32)> = (0..6).map(|i| (i, (i + 1) % 6)).collect();
    let mut tip = coords[0];
    for t in 0..4u32 {
        tip = [
            tip[0] + 1.25,
            tip[1],
            tip[2] + if t % 2 == 0 { 0.9 } else { -0.9 },
        ];
        coords.push(tip); // atoms 6, 7, 8, 9
        bonds.push(if t == 0 { (0, 6) } else { (5 + t, 6 + t) });
    }
    (coords, bonds)
}

/// Rebuilding a linear chain with the template's own torsion values must
/// reproduce the template coordinates exactly (spec Task 2 (i)).
#[test]
fn internal_roundtrip_linear() {
    let coords = zigzag_coords(12, 1.53);
    let tree = assert_roundtrip(&coords, &chain_bonds(12), "linear 12-chain");
    // Bonds 1-2 … 9-10 rotate; the two terminal bonds never do.
    assert_eq!(tree.n_vars(), 9, "linear 12-chain: 9 free torsions");
    assert_eq!(
        tree.n_steps(),
        9,
        "one step per free torsion, no rigid preamble"
    );
}

/// Same round-trip for a branched molecule: branch atoms are placed off the
/// backbone and must come back exactly, and the terminal branch bonds must
/// not add free variables.
#[test]
fn internal_roundtrip_branched() {
    let (coords, bonds) = branched_parts();
    let tree = assert_roundtrip(&coords, &bonds, "branched");
    // Backbone bonds 1-2 … 7-8 rotate; both chain-end bonds and both
    // one-atom branch bonds are terminal.
    assert_eq!(tree.n_vars(), 7, "branched: 7 free torsions");
}

/// Same round-trip for a molecule containing a ring: ring geometry (and the
/// ring-closure bond) must come back exactly, and cyclic bonds must not add
/// free variables.
#[test]
fn internal_roundtrip_ring() {
    let (coords, bonds) = ring_tail_parts();
    let tree = assert_roundtrip(&coords, &bonds, "ring+tail");
    // Only the tail bonds 0-6, 6-7 and 7-8 are free: the six ring bonds are
    // cyclic and the tail's last bond is terminal.
    assert_eq!(tree.n_vars(), 3, "ring+tail: 3 free torsions");
}

/// A molecule with no rotatable bonds (tetrahedral star: every bond is
/// terminal) still decomposes — zero variables, at least one step, exact
/// round-trip with an empty variable vector.
#[test]
fn internal_rigid_molecule_no_vars() {
    let s = 1.53 as F / (3.0 as F).sqrt();
    let coords = [[0.0, 0.0, 0.0], [s, s, s], [s, -s, -s], [-s, s, -s]];
    let bonds = [(0, 1), (0, 2), (0, 3)];
    let tree = assert_roundtrip(&coords, &bonds, "rigid star");
    assert_eq!(tree.n_vars(), 0, "terminal bonds are never rotatable");
    assert!(tree.n_steps() >= 1, "the rigid remainder is still a step");
}

/// THE ring-misclassification detector (spec ac-003(ii), Task 2 (ii)): under
/// ARBITRARY free-variable values, every template bond length (including the
/// ring-closure bond, which the round-trip test cannot see) and every bonded
/// angle must still match the template, and all coordinates must be finite.
/// Free torsions legitimately change and are not checked.
#[test]
fn internal_random_vars_preserve_bonded_geometry() {
    let pi = std::f64::consts::PI as F;
    let mut rng = SmallRng::seed_from_u64(20260828);
    for (what, (coords, bonds)) in [
        ("branched", branched_parts()),
        ("ring+tail", ring_tail_parts()),
    ] {
        let tree = InternalTree::from_frame(&frame_from_parts(&coords, &bonds))
            .unwrap_or_else(|e| panic!("{what}: template must decompose, got {e}"));
        assert!(
            tree.n_vars() > 0,
            "{what}: needs free torsions to randomize"
        );
        let angles = bonded_angle_triples(coords.len(), &bonds);
        for trial in 0..5 {
            let vars: Vec<F> = (0..tree.n_vars())
                .map(|_| (rng.random::<f64>() as F * 2.0 - 1.0) * pi)
                .collect();
            let rebuilt = rebuild_coords(&tree, &coords, &vars);
            for (i, p) in rebuilt.iter().enumerate() {
                assert!(
                    p.iter().all(|v| v.is_finite()),
                    "{what} trial {trial}: atom {i} is not finite: {p:?}"
                );
            }
            for &(i, j) in &bonds {
                let (i, j) = (i as usize, j as usize);
                let want = vdist(coords[i], coords[j]);
                let got = vdist(rebuilt[i], rebuilt[j]);
                assert!(
                    (got - want).abs() < 1e-9,
                    "{what} trial {trial}: bond ({i},{j}) = {got} vs template {want} — \
                     a ring-closure or otherwise non-free bond classified as a free \
                     variable breaks exactly this"
                );
            }
            for &(i, j, k) in &angles {
                let want = vangle(coords[i], coords[j], coords[k]);
                let got = vangle(rebuilt[i], rebuilt[j], rebuilt[k]);
                assert!(
                    (got - want).abs() < 1e-9,
                    "{what} trial {trial}: angle ({i},{j},{k}) = {got} rad vs template {want} rad"
                );
            }
        }
    }
}

/// The 1-`EXCLUDE_BONDS` (depth 3, AA convention) exclusion table: sorted,
/// self included, 1-2/1-3/1-4 partners in, 1-5 partners out.
#[test]
fn internal_exclusions_depth() {
    let tree =
        InternalTree::from_frame(&chain_frame(12, 1.53, true)).expect("linear chain decomposes");
    assert_eq!(
        tree.exclusions(0),
        &[0u32, 1, 2, 3][..],
        "atom 0 of a linear chain excludes exactly self + its 1-2/1-3/1-4 partners"
    );
    assert!(
        !tree.exclusions(0).contains(&4),
        "1-5 partners must be scored, not excluded — the chain must not thread itself"
    );
}

// ── Section: Task 3 — OverlapField: cell list over the periodic box ────────
//
// Review tests for `src/grow/field.rs` (spec Design §4c): insert/remove
// idempotence (retraction is O(1) removal, not a rebuild), minimum-image
// distances across the periodic boundary, 27-cell stencil de-duplication
// when the grid has fewer than 3 cells per axis, exclusion-list handling,
// and the `hard_scale` softening knob (spec §4f).

/// Insert / remove / re-insert must leave the field exactly consistent:
/// counts, occupancy, and what probe/nearest can see.
#[test]
fn field_insert_remove_idempotent() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0; 10],
        (0..10).collect(),
        vec![0; 10],
        2.5,
    );

    let placed: [[F; 3]; 5] = [
        [1.0, 1.0, 1.0],
        [4.0, 4.0, 4.0],
        [7.0, 7.0, 7.0],
        [1.0, 7.0, 1.0],
        [7.0, 1.0, 4.0],
    ];
    for (s, &p) in placed.iter().enumerate() {
        field.insert(s, p);
    }
    assert_eq!(field.n_placed(), 5);
    for (s, &p) in placed.iter().enumerate() {
        assert!(field.is_placed(s), "slot {s} must be placed");
        assert_eq!(field.position(s), p, "slot {s} position");
    }
    for s in 5..10 {
        assert!(!field.is_placed(s), "slot {s} was never inserted");
    }

    field.remove(1);
    field.remove(3);
    assert_eq!(field.n_placed(), 3);
    assert!(!field.is_placed(1) && !field.is_placed(3));
    field.remove(1); // documented no-op on an empty slot
    assert_eq!(field.n_placed(), 3, "double remove must be a no-op");

    // Re-insert the two slots elsewhere: probe/nearest must see the new
    // positions and have forgotten the old ones.
    field.insert(1, [2.2, 1.0, 1.0]);
    field.insert(3, [9.4, 1.0, 1.0]);
    assert_eq!(field.n_placed(), 5);
    assert!(field.is_placed(1) && field.is_placed(3));
    let d = field.nearest(9, [2.2, 1.0, 1.0], &[]);
    assert!(
        d < 1e-12,
        "nearest at slot 1's new position = {d}, must be 0"
    );
    assert_eq!(
        field.probe(9, [4.0, 4.0, 4.0], &[], 1.0, 0.5, 1e9),
        Probe::Room(0.0),
        "slot 1's old position must be empty after re-insertion"
    );
    assert_eq!(
        field.probe(9, [2.2, 1.0, 1.0], &[], 1.0, 0.5, 1e9),
        Probe::Blocked,
        "slot 1's new position must be hard-blocked"
    );

    for s in 0..5 {
        field.remove(s);
    }
    assert_eq!(field.n_placed(), 0);
    for p in [[0.0; 3], [5.0; 3], [9.9, 0.1, 5.0]] {
        assert_eq!(
            field.probe(9, p, &[], 1.0, 0.5, 1e9),
            Probe::Room(0.0),
            "an empty field must answer Room(0.0) everywhere"
        );
    }
}

/// Distances must use the minimum image: 0.5 and 9.9 in a 10 Å periodic box
/// are 0.6 Å apart, well inside the 2.0 Å contact.
#[test]
fn field_minimum_image_across_boundary() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0, 1.0],
        vec![0, 1],
        vec![0, 0],
        2.5,
    );
    field.insert(0, [0.5, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [9.9, 5.0, 5.0], &[], 1.0, 0.4, 1e9),
        Probe::Blocked,
        "0.6 Å across the periodic boundary is inside the 2.0 Å contact"
    );
    let d = field.nearest(1, [9.9, 5.0, 5.0], &[]);
    assert!(
        (d - 0.6).abs() < 1e-12,
        "nearest across the boundary = {d}, expected 0.6"
    );
}

/// Box 4³ with cutoff 2.5 → a single cell per axis: all 27 stencil offsets
/// alias the same cell, and a naive stencil would count the one neighbour up
/// to 27 times. The penalty is hand-derived and binary-exact, so any double
/// counting shows up as an integer multiple.
#[test]
fn field_small_box_stencil_dedup() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [4.0; 3],
        [true; 3],
        vec![0.5, 0.5],
        vec![0, 1],
        vec![0, 0],
        2.5,
    );
    field.insert(0, [1.0, 2.0, 2.0]);
    // d = 1.5, contact = 1.0, shell = contact + soft_shell = 2.0. Exactly one
    // pair inside the soft shell: gap = (shell² − d²)/shell² = (4 − 2.25)/4
    // = 0.4375, penalty = gap² = 0.19140625.
    let probe = field.probe(1, [2.5, 2.0, 2.0], &[], 1.0, 1.0, 1e9);
    let Probe::Room(penalty) = probe else {
        panic!("d = 1.5 > contact = 1.0 must not be Blocked, got {probe:?}");
    };
    assert!(
        (penalty - 0.19140625).abs() < 1e-12,
        "penalty = {penalty}, expected exactly 0.19140625 — 2× means the \
         aliased cell was visited twice"
    );
    let d = field.nearest(1, [2.5, 2.0, 2.0], &[]);
    assert!((d - 1.5).abs() < 1e-12, "nearest = {d}, expected 1.5");
}

/// Same-molecule template atoms on the exclusion list are neither blocking
/// nor charged; a same-molecule atom OFF the list still blocks.
#[test]
fn field_exclusions_respected() {
    // Three slots of the same molecule (mol 7): template atoms 0, 1 and 4.
    let mut field = OverlapField::new(
        [0.0; 3],
        [20.0; 3],
        [true; 3],
        vec![1.0; 3],
        vec![7; 3],
        vec![0, 1, 4],
        2.5,
    );
    let excluded: &[u32] = &[0, 1, 2, 3]; // template atom 4 is NOT excluded
    field.insert(0, [5.0, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [5.5, 5.0, 5.0], excluded, 1.0, 0.5, 1e9),
        Probe::Room(0.0),
        "deep inside slot 0's hard core, but its template atom is excluded → Room"
    );
    field.insert(2, [6.5, 5.0, 5.0]);
    assert_eq!(
        field.probe(1, [6.0, 5.0, 5.0], excluded, 1.0, 0.5, 1e9),
        Probe::Blocked,
        "template atom 4 is not on the exclusion list and violates → Blocked"
    );
}

/// `hard_scale` shrinks the refusal core (the §4f soften fallback): a pair at
/// d = 1.8 with contact 2.0 is Blocked at scale 1.0 and merely charged at 0.8.
#[test]
fn field_hard_scale_softens() {
    let mut field = OverlapField::new(
        [0.0; 3],
        [10.0; 3],
        [true; 3],
        vec![1.0, 1.0],
        vec![0, 1],
        vec![0, 0],
        3.0,
    );
    field.insert(0, [5.0, 5.0, 5.0]);
    let p = [6.8, 5.0, 5.0];
    assert_eq!(
        field.probe(1, p, &[], 1.0, 0.4, 1e9),
        Probe::Blocked,
        "hard_scale 1.0: d = 1.8 < contact = 2.0"
    );
    match field.probe(1, p, &[], 0.8, 0.4, 1e9) {
        Probe::Room(penalty) => assert!(
            penalty > 0.0,
            "d = 1.8 sits inside the softened soft shell (1.6..2.0): allowed \
             but charged, got penalty {penalty}"
        ),
        Probe::Blocked => panic!("hard_scale 0.8 must admit d = 1.8 > 1.6"),
    }
}

// ── Section: Task 4 — TorsionPrior: sampling + C∞ calibration ──────────────
//
// RED via compile failure until `src/grow/prior.rs` grows the frozen Task 4
// surface (spec Design §4a′):
//
//   TorsionPrior::sample(&self, template_value: F, rng: &mut impl rand::Rng) -> F
//   TorsionPrior::three_state_from_c_inf(c_inf: F, theta: F) -> TorsionPrior
//
// `sample` semantics: Uniform → uniform on (−π, π]; Template{kappa} →
// von-Mises-like around `template_value`; States → pick a state angle by
// normalized weight, exact angle, no jitter (v1).

/// Mean characteristic ratio C_n = ⟨R²⟩ / ((n−1)·b²) over `n_samples`
/// conformers of an isolated chain (end-to-end between the graph ends,
/// no PBC), every free torsion drawn independently from `prior`.
fn sampled_c_n(prior: &TorsionPrior, n_beads: usize, bond: F, n_samples: usize, seed: u64) -> F {
    let template = zigzag_coords(n_beads, bond);
    let tree = InternalTree::from_frame(&frame_from_parts(&template, &chain_bonds(n_beads)))
        .expect("chain template decomposes");
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut vars = vec![0.0 as F; tree.n_vars()];
    let mut sum_r2 = 0.0 as F;
    for _ in 0..n_samples {
        for k in 0..tree.n_steps() {
            if let Some(v) = tree.step_var(k) {
                vars[v] = prior.sample(tree.template_var(k), &mut rng);
            }
        }
        let rebuilt = rebuild_coords(&tree, &template, &vars);
        let r = vsub(rebuilt[n_beads - 1], rebuilt[0]);
        sum_r2 += vdot(r, r);
    }
    sum_r2 / n_samples as F / ((n_beads - 1) as F * bond * bond)
}

/// `three_state_from_c_inf(5.5, 109.47°)` — the PEO calibration (spec §5.1):
/// States with trans at |φ| = π and gauche± at ±π/3 (IUPAC absolute
/// convention), weights (p_t, p_g, p_g) normalized, p_t ≈ 0.645 solved from
/// C∞ = C_FRC·(1+⟨cosφ'⟩)/(1−⟨cosφ'⟩) with φ' measured from trans.
#[test]
fn prior_three_state_calibration() {
    let pi = std::f64::consts::PI as F;
    let theta = (109.47 as F).to_radians();
    let prior = TorsionPrior::three_state_from_c_inf(5.5, theta);
    let TorsionPrior::States(states) = &prior else {
        panic!("three_state_from_c_inf must return TorsionPrior::States, got {prior:?}");
    };
    assert_eq!(states.len(), 3, "three states: trans + gauche±");
    let total: F = states.iter().map(|&(_, w)| w).sum();
    assert!(
        (total - 1.0).abs() < 1e-12,
        "weights must sum to 1, got {total}"
    );

    let (trans, mut gauche): (Vec<_>, Vec<_>) = states
        .iter()
        .copied()
        .partition(|&(a, _)| (a.abs() - pi).abs() < 1e-9);
    assert_eq!(
        trans.len(),
        1,
        "exactly one trans state at |φ| = π, got {states:?}"
    );
    // C_FRC = (1−cosθ)/(1+cosθ) ≈ 2.000; r = 5.5/C_FRC ≈ 2.750;
    // x = (r−1)/(r+1) ≈ 0.4667; p_t = (2x+1)/3 ≈ 0.6445.
    assert!(
        (trans[0].1 - 0.645).abs() < 0.01,
        "trans weight = {}, expected ≈ 0.645 for C∞ = 5.5",
        trans[0].1
    );
    gauche.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("finite angles"));
    assert_eq!(gauche.len(), 2, "two gauche states");
    assert!(
        (gauche[0].0 + pi / 3.0).abs() < 1e-9,
        "gauche− at −π/3, got {}",
        gauche[0].0
    );
    assert!(
        (gauche[1].0 - pi / 3.0).abs() < 1e-9,
        "gauche+ at +π/3, got {}",
        gauche[1].0
    );
    assert!(
        (gauche[0].1 - gauche[1].1).abs() < 1e-12,
        "gauche weights must be equal: {} vs {}",
        gauche[0].1,
        gauche[1].1
    );
}

/// The ac-006(1) regression baseline: uniform torsions on a fixed-109.5°
/// chain are the freely rotating chain, C∞ = (1−cosθ)/(1+cosθ) = 2.00
/// exactly. This pins the NeRF rebuild + sampler combination to an analytic
/// value; the fixed seed makes the sampled number stable.
#[test]
fn prior_uniform_freely_rotating_c_inf() {
    let c_n = sampled_c_n(&TorsionPrior::Uniform, 200, 1.53, 600, 42);
    assert!(
        (c_n - 2.0).abs() <= 0.1,
        "freely rotating chain: C_n = {c_n:.4}, analytic C∞ = 2.00 ± 0.1 (spec §5.1)"
    );
}

/// `States([(π, 1)])` is the all-trans prior, and the zigzag template IS the
/// all-trans conformer — sampling must reproduce the template coordinates.
#[test]
fn prior_states_trans_only_is_all_trans() {
    let pi = std::f64::consts::PI as F;
    let template = zigzag_coords(20, 1.53);
    let tree = InternalTree::from_frame(&frame_from_parts(&template, &chain_bonds(20)))
        .expect("chain template decomposes");
    let prior = TorsionPrior::States(vec![(pi, 1.0)]);
    let mut rng = SmallRng::seed_from_u64(7);
    let mut vars = vec![0.0 as F; tree.n_vars()];
    for k in 0..tree.n_steps() {
        if let Some(v) = tree.step_var(k) {
            vars[v] = prior.sample(tree.template_var(k), &mut rng);
        }
    }
    let rebuilt = rebuild_coords(&tree, &template, &vars);
    let dev = max_abs_dev(&template, &rebuilt);
    assert!(
        dev < 1e-6,
        "an all-trans States prior must reproduce the (already all-trans) \
         zigzag template: ‖Δ‖∞ = {dev:e}"
    );
}

/// The calibration round-trip (spec chain-statistics ladder, step 3): a chain
/// sampled from `three_state_from_c_inf(5.5, 109.47°)` must measure back
/// C_n = 5.5 ± 0.3. Finite-n droop is real at n = 200; if this fails
/// marginally when Task 4 lands, raise n_beads to 400 and keep the tolerance.
#[test]
fn prior_ris_calibrated_c_inf() {
    let theta = (109.47 as F).to_radians();
    let prior = TorsionPrior::three_state_from_c_inf(5.5, theta);
    let c_n = sampled_c_n(&prior, 200, 1.53, 600, 42);
    assert!(
        (c_n - 5.5).abs() <= 0.3,
        "RIS-calibrated chain: C_n = {c_n:.4}, target C∞ = 5.5 ± 0.3"
    );
}

// ── Section: Task 5 — GrowStage: constructive growth end-to-end ───────────
//
// Exercises `src/grow/driver.rs` (spec Design §4, Task 5): the `GrowConfig`
// builder surface, the `GrowError::{NoBox, TriclinicCell, FixedTarget}`
// rejections, and the constructive all-grow pack itself.
//
// Mixed grow+gencan composition in one call is gone with the monolithic
// entry (engine-entry-split); the explicit chain over a fixed matrix is
// covered by `grow_then_gencan_chaining_over_fixed_matrix` below.

/// Minimum-image displacement in a cubic periodic box of edge `l`.
fn min_image_cubic(d: [F; 3], l: F) -> [F; 3] {
    std::array::from_fn(|k| d[k] - (d[k] / l).round() * l)
}

/// Radius of gyration of one linear chain, bond-graph-unwrapped under the
/// minimum image so a chain wrapping the periodic box is not inflated
/// (`examples/pack_peo/main.rs::radius_of_gyration` pattern). For a linear
/// `(i, i+1)` bond list the BFS unwrap IS the sequential walk below.
fn linear_chain_rg(xyz: &[[F; 3]], l: F) -> F {
    let mut un: Vec<[F; 3]> = Vec::with_capacity(xyz.len());
    un.push(xyz[0]);
    for w in xyz.windows(2) {
        let prev = *un.last().expect("unwrap walk is never empty");
        let d = min_image_cubic(vsub(w[1], w[0]), l);
        un.push([prev[0] + d[0], prev[1] + d[1], prev[2] + d[2]]);
    }
    let n = un.len() as F;
    let mut com = [0.0 as F; 3];
    for p in &un {
        for k in 0..3 {
            com[k] += p[k] / n;
        }
    }
    let msd: F = un
        .iter()
        .map(|&p| {
            let d = vsub(p, com);
            vdot(d, d)
        })
        .sum::<F>()
        / n;
    msd.sqrt()
}

/// Brute-force O(N²) minimum INTER-molecular atom distance under the minimum
/// image. `mol_of[i]` is atom i's molecule id; same-molecule pairs are the
/// exclusion table's jurisdiction and are skipped.
fn min_inter_distance(pos: &[[F; 3]], mol_of: &[usize], l: F) -> F {
    let mut min = F::INFINITY;
    for i in 0..pos.len() {
        for j in i + 1..pos.len() {
            if mol_of[i] == mol_of[j] {
                continue;
            }
            let d = min_image_cubic(vsub(pos[j], pos[i]), l);
            min = min.min(vdot(d, d).sqrt());
        }
    }
    min
}

/// Number of INTER-molecular atom pairs closer than `d` under the minimum
/// image — the companion counter to [`min_inter_distance`]: the minimum
/// answers "how bad is the worst pair", this answers "how many pairs are
/// bad at all", so a single lucky-looking minimum cannot hide a population
/// of violations. Cell-binned (edge `d`, 27-cell stencil) so it stays cheap
/// on a 15 000-atom melt; the O(N²) walk is kept as the fallback whenever
/// the box holds fewer than 3 cells per axis, where the stencil would alias
/// a cell onto itself.
fn count_inter_pairs_below(pos: &[[F; 3]], mol_of: &[usize], l: F, d: F) -> usize {
    let nc = (l / d).floor() as isize;
    if nc < 3 {
        let mut n = 0usize;
        for i in 0..pos.len() {
            for j in i + 1..pos.len() {
                if mol_of[i] == mol_of[j] {
                    continue;
                }
                let s = min_image_cubic(vsub(pos[j], pos[i]), l);
                if vdot(s, s).sqrt() < d {
                    n += 1;
                }
            }
        }
        return n;
    }

    let edge = l / nc as F;
    let cell_of = |p: [F; 3]| -> [isize; 3] {
        std::array::from_fn(|k| {
            let w = p[k] - (p[k] / l).floor() * l;
            ((w / edge) as isize).clamp(0, nc - 1)
        })
    };
    let index = |c: [isize; 3]| -> usize {
        ((c[0].rem_euclid(nc) * nc + c[1].rem_euclid(nc)) * nc + c[2].rem_euclid(nc)) as usize
    };

    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); (nc * nc * nc) as usize];
    for (i, &p) in pos.iter().enumerate() {
        bins[index(cell_of(p))].push(i);
    }

    let mut n = 0usize;
    for (i, &p) in pos.iter().enumerate() {
        let c = cell_of(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    for &j in &bins[index([c[0] + dx, c[1] + dy, c[2] + dz])] {
                        if j <= i || mol_of[i] == mol_of[j] {
                            continue;
                        }
                        let s = min_image_cubic(vsub(pos[j], p), l);
                        if vdot(s, s).sqrt() < d {
                            n += 1;
                        }
                    }
                }
            }
        }
    }
    n
}

/// One all-`Grow` pack: `copies` copies of an `n_beads` zigzag chain
/// (bond 1.53, bonds present, `TorsionPrior::Uniform`) in a cubic periodic
/// box `[0, box_len]³`, tolerance 2.0, 50 loops.
fn grow_pack(
    seed: u64,
    copies: usize,
    n_beads: usize,
    box_len: F,
) -> Result<PackResult, PackError> {
    let target = Target::new(chain_frame(n_beads, 1.53, true), copies);
    CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [box_len; 3], [true; 3])
        .run(&[target], 50)
}

/// The Task 5 headline: 8 × 12-bead chains at moderate fill grow to a
/// feasible state CONSTRUCTIVELY — `fdist` is strict bitwise zero because
/// hard-core violation is rejected, never penalized (spec §4c, ac-004) — and
/// the copies are distinct conformers (ac-006's anti-pathology).
#[test]
fn grow_packs_small_melt_constructively() {
    let (copies, n_beads, l) = (8usize, 12usize, 26.0 as F);
    let result = grow_pack(7, copies, n_beads, l)
        .expect("an all-Grow pack must run the growth solver and succeed (Task 5)");

    assert_eq!(
        result.fdist.to_bits(),
        (0.0 as F).to_bits(),
        "fdist = {} — hard rejection makes fdist == 0.0 a constructive \
         guarantee, not a convergence hope (ac-004)",
        result.fdist
    );
    assert_eq!(
        result.softened, 0,
        "moderate fill must not need the softening fallback"
    );
    assert!(result.converged, "softened == 0 growth must be converged");

    let pos = result.positions();
    assert_eq!(pos.len(), copies * n_beads, "8 copies × 12 beads");

    // Independent ruler: brute-force inter-molecular minimum distance.
    let mol_of: Vec<usize> = (0..pos.len()).map(|i| i / n_beads).collect();
    let dmin = min_inter_distance(&pos, &mol_of, l);
    assert!(
        dmin >= 2.0 - 1e-9,
        "minimum inter-molecular atom distance = {dmin}, must be ≥ tolerance \
         2.0 − 1e-9 under the minimum image"
    );

    // The copies must be distinct conformers, not 8 clones of one shape.
    let rgs: Vec<F> = (0..copies)
        .map(|c| linear_chain_rg(&pos[c * n_beads..(c + 1) * n_beads], l))
        .collect();
    let (lo, hi) = rgs
        .iter()
        .fold((F::INFINITY, F::NEG_INFINITY), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
    assert!(
        hi - lo > 1e-6,
        "per-copy Rg all within [{lo}, {hi}] — copies must be DISTINCT \
         conformers (ac-006 anti-pathology)"
    );
}

/// Same seed, same inputs → bitwise-identical positions (spec determinism
/// gate).
#[test]
fn grow_deterministic_same_seed() {
    let a = grow_pack(7, 8, 12, 26.0).expect("run A must succeed");
    let b = grow_pack(7, 8, 12, 26.0).expect("run B must succeed");
    let (pa, pb) = (a.positions(), b.positions());
    assert_eq!(pa.len(), pb.len(), "same atom count");
    for (i, (p, q)) in pa.iter().zip(&pb).enumerate() {
        for k in 0..3 {
            assert_eq!(
                p[k].to_bits(),
                q[k].to_bits(),
                "atom {i} axis {k}: same seed must reproduce bitwise \
                 ({} vs {})",
                p[k],
                q[k]
            );
        }
    }
}

/// Different seeds must produce different packs — a cheap sanity check that
/// the seed is actually consumed by the growth RNG streams.
#[test]
fn grow_seed_changes_result() {
    let a = grow_pack(7, 8, 12, 26.0).expect("seed 7 must succeed");
    let b = grow_pack(8, 8, 12, 26.0).expect("seed 8 must succeed");
    let (pa, pb) = (a.positions(), b.positions());
    assert_eq!(pa.len(), pb.len(), "same atom count");
    let differs = pa
        .iter()
        .zip(&pb)
        .any(|(p, q)| (0..3).any(|k| p[k].to_bits() != q[k].to_bits()));
    assert!(
        differs,
        "seeds 7 and 8 produced bitwise-identical packs — the seed is not \
         being consumed"
    );
}

/// THE ac-012 detector. Box edge 200 Å: two 12-bead chains can never
/// interact, so copy 0's growth must not depend on copy 1's existence.
/// Per-(copy, step) hashed RNG streams + round-snapshot semantics (spec §2)
/// make copy 0's atoms bitwise identical between a 1-copy and a 2-copy run;
/// a single global RNG stream fails exactly this.
#[test]
fn grow_copy_stream_independence() {
    let a = grow_pack(7, 1, 12, 200.0).expect("1-copy run must succeed");
    let b = grow_pack(7, 2, 12, 200.0).expect("2-copy run must succeed");
    let (pa, pb) = (a.positions(), b.positions());
    assert_eq!(pa.len(), 12, "run A: 1 copy × 12 beads");
    assert_eq!(pb.len(), 24, "run B: 2 copies × 12 beads");
    for i in 0..12 {
        for k in 0..3 {
            assert_eq!(
                pa[i][k].to_bits(),
                pb[i][k].to_bits(),
                "copy 0 atom {i} axis {k}: {} vs {} — copy 0's stream must \
                 be independent of copy 1's existence (hashed per-(copy, \
                 step) RNG, spec §2 / ac-012)",
                pa[i][k],
                pb[i][k]
            );
        }
    }
}

/// Two different Grow species in one pack: both grow, both constructive,
/// both present in target order.
#[test]
fn grow_two_species() {
    let long = Target::new(chain_frame(12, 1.53, true), 4);
    let short = Target::new(chain_frame(6, 1.53, true), 6);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .run(&[long, short], 50)
        .expect("an all-Grow two-species pack must succeed");
    assert_eq!(
        result.fdist.to_bits(),
        (0.0 as F).to_bits(),
        "constructive fdist == 0.0 across species, got {}",
        result.fdist
    );
    assert_eq!(result.softened, 0, "no softening at moderate fill");
    assert_eq!(
        result.positions().len(),
        4 * 12 + 6 * 6,
        "both species' copies must be present"
    );
}

/// A Grow target with neither a periodic box nor a cell is a named `NoBox`
/// error — growth needs a final volume from atom 0 (spec §4b/§7), and
/// molpack does not guess one.
#[test]
fn grow_requires_box() {
    let target = Target::new(chain_frame(12, 1.53, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .run(&[target], 50)
        .expect_err("Grow without a periodic box or cell must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::NoBox,
                ..
            }
        ),
        "expected PackError::Grow(NoBox), got {err:?}"
    );
}

/// A Grow target in a non-orthorhombic cell is a named `TriclinicCell` error
/// — the v1 `OverlapField` is orthorhombic-only (spec §10, Out of scope).
#[test]
fn grow_rejects_triclinic() {
    let target = Target::new(chain_frame(12, 1.53, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_cell([20.0, 20.0, 20.0], [60.0, 60.0, 60.0], [true; 3])
        .run(&[target], 50)
        .expect_err("Grow in a 60°/60°/60° cell must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::TriclinicCell,
                ..
            }
        ),
        "expected PackError::Grow(TriclinicCell), got {err:?}"
    );
}

/// `fixed_at` combined with Grow is a named `FixedTarget` error. The box is
/// valid, so the fixed placement is the only reason to refuse.
#[test]
fn grow_rejects_fixed_target() {
    let target = Target::new(chain_frame(12, 1.53, true), 1).fixed_at([5.0, 5.0, 5.0]);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .run(&[target], 50)
        .expect_err("fixed_at combined with Grow must be rejected");
    assert!(
        matches!(
            err,
            PackError::Grow {
                source: GrowError::FixedTarget,
                ..
            }
        ),
        "expected PackError::Grow(FixedTarget), got {err:?}"
    );
}

/// Compile-level pin of the frozen `GrowConfig` builder surface (consuming
/// builders, `Target`-style; every knob has a default — values are NOT
/// asserted here, only the method names and signatures).
#[test]
fn grow_config_builder_chain() {
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_trials(16)
        .with_selectivity(1.5)
        .with_soft_shell(0.8)
        .with_retract(6)
        .with_relax(20, 4)
        .with_soften_after(30)
        .with_min_hard_scale(0.85)
        .with_exclusion_depth(2);
    let dbg = format!("{cfg:?}");
    assert!(!dbg.is_empty(), "GrowConfig must keep a useful Debug impl");
}

// ── Section: Tasks 6-7 — restraint hard rejection + handler wiring ─────────
//
// runtime-RED until Task 6 (growth currently IGNORES `AtomRestraint`s during
// placement, so a restrained Grow pack scatters atoms outside the region and
// reports frest > 0) and Task 7 (growth currently emits NO `on_step` events
// and never polls `should_stop`).
//
// Task 6 contract (spec Design §3, ac-007): every candidate atom position is
// checked against the target's restraints via the existing
// `AtomRestraint::f`; `f > 0` is a hard rejection, same treatment as a
// hard-core violation. `frest == 0.0` thereby becomes a CONSTRUCTIVE
// guarantee, exactly like `fdist == 0.0` — strict zero, not `< precision`.
//
// Task 7 contract (spec Design §4g): one `StepInfo` per growth round with
// `loop_idx` = round number (1-based, strictly increasing), `radscale` =
// current hard-core scale (1.0 while unsoftened), and fdist/frest = 0.0
// while the hard-rejection regime holds. `Handler::should_stop() == true`
// aborts growth: `pack` still returns Ok, with `converged == false`.

/// One recorded `on_step` event: `(loop_idx, radscale, fdist, frest)`,
/// shared with the test through an `Arc` so it stays readable after `pack`
/// consumes the handler.
type StepEvents = Arc<Mutex<Vec<(usize, F, F, F)>>>;

/// Records every growth `on_step` event and (for the abort test) requests a
/// stop once `stop_after` events have been seen (`usize::MAX` = never stop).
struct Recorder {
    events: StepEvents,
    stop_after: usize,
}

impl Handler for Recorder {
    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        self.events.lock().expect("recorder mutex").push((
            info.loop_idx,
            info.radscale,
            info.fdist,
            info.frest,
        ));
    }

    fn should_stop(&self) -> bool {
        self.events.lock().expect("recorder mutex").len() >= self.stop_after
    }
}

/// The softening contract, read off a recorded `radscale` trajectory.
///
/// `radscale` on a growth `StepInfo` IS the driver's hard-core scale, so the
/// recorded rounds are the only public window on the softening schedule. The
/// schedule is a **ladder**: monotone non-increasing, bounded below by the
/// configured floor, and dropping by at most ONE 3 % rung per round. One rung
/// per round is not a style preference — a rung is earned by some chain
/// accumulating `soften_after` dead ends, and a chain can dead-end at most
/// once per round (`src/grow/driver.rs` gives every pending chain exactly one
/// commit attempt per round), so two rungs inside one round prove that some
/// path shrank the core without going through the per-chain counter. That is
/// the free-fall of debt D-01 (`.claude/notes/notes.md`, 2026-09-02).
fn assert_radscale_ladder(events: &[(usize, F, F, F)], floor: F, what: &str) {
    assert!(
        !events.is_empty(),
        "{what}: growth emitted no on_step events, so the softening schedule \
         is unobservable — every growth round must emit one StepInfo"
    );
    for &(round, radscale, _, _) in events {
        assert!(
            radscale <= 1.0 + 1e-12 && radscale >= floor - 1e-12,
            "{what}: round {round} reports radscale = {radscale}, outside \
             [{floor}, 1.0] — the hard-core scale starts at 1.0 and the \
             softening floor (min_hard_scale × 1) is a hard bound, not a hint"
        );
    }
    for w in events.windows(2) {
        let ((prev_round, prev), (round, cur)) = ((w[0].0, w[0].1), (w[1].0, w[1].1));
        assert!(
            cur <= prev + 1e-12,
            "{what}: radscale rose {prev} → {cur} between rounds {prev_round} \
             and {round} — softening is monotone (recovery is grow-axes work, \
             not something the driver may do silently)"
        );
        assert!(
            cur >= prev * 0.97 - 1e-12,
            "{what}: radscale fell {prev} → {cur} between rounds {prev_round} \
             and {round}, more than ONE 3 % rung in a single round. A chain \
             dead-ends at most once per round, so every rung must be earned by \
             soften_after dead ends on some chain; a multi-rung round is the \
             per-dead-end free-fall of debt D-01, which the now-removed \
             regrow_budget path in driver.rs used to produce"
        );
    }
}

/// THE ac-007 test: growth under an `InsideSphereRestraint` hard-rejects any
/// candidate with `f > 0`, so EVERY atom lands inside the sphere and
/// `frest == 0.0` is a strict constructive zero — same status as `fdist`.
/// The radius-10 sphere at the box center is strictly interior to [0, 30]³,
/// so the raw (unwrapped) distance is the right ruler.
#[test]
fn grow_restraint_hard_rejects() {
    let (center, radius) = ([15.0 as F, 15.0, 15.0], 10.0 as F);
    let target = Target::new(chain_frame(6, 1.53, true), 4)
        .with_restraint(InsideSphereRestraint::new(center, radius));
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
        .run(&[target], 50)
        .expect("a Grow pack with a generously feasible sphere restraint must succeed");

    let pos = result.positions();
    assert_eq!(pos.len(), 4 * 6, "4 copies × 6 beads");
    for (i, &p) in pos.iter().enumerate() {
        let d = vdist(p, center);
        assert!(
            d <= radius + 1e-9,
            "atom {i} at distance {d} from the sphere center — candidates \
             violating a restraint must be hard-rejected during growth \
             (spec §3, ac-007), never placed and penalized later"
        );
    }
    assert_eq!(
        result.frest.to_bits(),
        (0.0 as F).to_bits(),
        "frest = {} — hard rejection makes frest == 0.0 a constructive \
         guarantee, strict zero, not < precision (ac-007)",
        result.frest
    );
    assert_eq!(
        result.fdist.to_bits(),
        (0.0 as F).to_bits(),
        "fdist = {} — the hard-core guarantee must survive restraint wiring",
        result.fdist
    );
    assert_eq!(
        result.softened, 0,
        "24 beads in a radius-10 sphere is generous — no softening fallback"
    );
}

/// An unsatisfiable restraint must not hang and must not be reported as a
/// success: a radius-1.0 sphere cannot hold a 6-bead chain (the 1-3 distance
/// alone is ~2.5 Å > the sphere's 2.0 Å diameter, at any torsion). The
/// driver's forced-placement escape must surface through the result:
/// Ok, `converged == false`, `softened > 0`.
#[test]
fn grow_restraint_infeasible_is_not_silent() {
    let target = Target::new(chain_frame(6, 1.53, true), 2)
        .with_restraint(InsideSphereRestraint::new([15.0, 15.0, 15.0], 1.0));
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
        .run(&[target], 3)
        .expect(
            "an infeasible restraint is not an Err — the budgeted escape \
             reports through converged/softened",
        );
    assert!(
        !result.converged,
        "a pack that could not satisfy its restraint must not report converged"
    );
    assert!(
        result.softened > 0,
        "softened = 0 on an unsatisfiable restraint — the constructive \
         guarantee was silently broken instead of counted (spec §4f)"
    );
}

/// Task 7 headline: the small-melt setup of
/// `grow_packs_small_melt_constructively` must emit one `on_step` per growth
/// round — 1-based strictly increasing `loop_idx`, `radscale == 1.0`
/// throughout (no softening in this easy box), and fdist/frest == 0.0 while
/// hard rejection holds (spec §4g).
#[test]
fn grow_emits_step_events() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let target = Target::new(chain_frame(12, 1.53, true), 8);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: usize::MAX,
        }))
        .run(&[target], 50)
        .expect("the Task 5 small-melt baseline must still succeed");
    // Guards the radscale claim below: this box needs no softening.
    assert_eq!(result.softened, 0, "easy box must not soften");

    let events = events.lock().expect("recorder mutex");
    assert!(
        !events.is_empty(),
        "growth emitted no on_step events — every growth round must emit one \
         StepInfo (spec §4g, Task 7)"
    );
    assert_eq!(
        events[0].0, 1,
        "growth rounds are 1-based: first loop_idx must be 1, got {}",
        events[0].0
    );
    for w in events.windows(2) {
        assert!(
            w[1].0 > w[0].0,
            "loop_idx must be strictly increasing: {:?} then {:?}",
            w[0],
            w[1]
        );
    }
    for &(round, radscale, fdist, frest) in events.iter() {
        assert_eq!(
            radscale, 1.0,
            "round {round}: unsoftened growth must report radscale == 1.0 \
             (radscale = current hard-core scale, spec §4g), got {radscale}"
        );
        assert_eq!(
            fdist.to_bits(),
            (0.0 as F).to_bits(),
            "round {round}: fdist = {fdist} — hard rejection holds, so \
             per-round fdist is 0.0 by construction (spec §4g)"
        );
        assert_eq!(
            frest.to_bits(),
            (0.0 as F).to_bits(),
            "round {round}: frest = {frest} — hard rejection holds, so \
             per-round frest is 0.0 by construction (spec §4g)"
        );
    }
}

/// `Handler::should_stop() == true` (here: after the 2nd on_step) aborts
/// growth: `pack` still returns Ok, `converged == false`, and only a handful
/// of rounds ran instead of the full schedule.
#[test]
fn grow_should_stop_aborts() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let target = Target::new(chain_frame(12, 1.53, true), 8);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: 2,
        }))
        .run(&[target], 50)
        .expect("a should_stop abort still returns Ok");
    assert!(
        !result.converged,
        "an aborted grow pack must not be reported as converged"
    );
    let n = events.lock().expect("recorder mutex").len();
    assert!(
        n >= 2,
        "growth emitted only {n} on_step events — should_stop can only \
         trigger if the round events of Task 7 actually fire"
    );
    assert!(
        n < 10,
        "growth recorded {n} rounds after the stop request — should_stop \
         must abort the growth loop, not be ignored until completion"
    );
}

/// Abort used to write unplaced atoms as `[0,0,0]` then stamp the full
/// template topology on top: 0-length bonds (adjacent sentinels) and
/// box-scale bonds (placed atom ↔ origin). The assembled 1-2 distances
/// must stay the template bond even when `converged == false`.
#[test]
fn grow_abort_keeps_bonded_geometry() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let (n_beads, copies, bond, box_len) = (12usize, 8usize, 1.53 as F, 26.0 as F);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [box_len; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: 2,
        }))
        .run(&[target], 50)
        .expect("a should_stop abort still returns Ok");
    assert!(
        !result.converged,
        "an aborted grow pack must not be reported as converged"
    );
    let n = events.lock().expect("recorder mutex").len();
    assert!(
        n < 10,
        "growth recorded {n} rounds — this test is only meaningful on an abort"
    );

    let pos = result.positions();
    assert_eq!(pos.len(), copies * n_beads);
    let bonds = result
        .frame
        .get("bonds")
        .expect("assembled frame keeps bonds");
    let ai = bonds.get_uint("atomi").expect("atomi");
    let aj = bonds.get_uint("atomj").expect("atomj");
    assert!(
        !ai.is_empty(),
        "template 1-2 bonds must be tiled onto copies"
    );

    let mut lo = F::INFINITY;
    let mut hi = F::NEG_INFINITY;
    for (&a, &b) in ai.iter().zip(aj.iter()) {
        let d = min_image_cubic(vsub(pos[a as usize], pos[b as usize]), box_len);
        let len = vdot(d, d).sqrt();
        lo = lo.min(len);
        hi = hi.max(len);
    }
    assert!(
        lo > 0.5,
        "bonded min {lo} Å — unplaced origin sentinels produce 0-length bonds"
    );
    assert!(
        hi < 5.0,
        "bonded max {hi} Å — a placed atom bonded to an origin sentinel \
         is box-scale, not a chemical 1-2 distance"
    );
    assert!(
        (lo - bond).abs() < 0.05 && (hi - bond).abs() < 0.05,
        "bonded range [{lo}, {hi}] Å must stay the template bond {bond} Å"
    );
}

/// Characterization golden for the abort completion path
/// (`.claude/specs/stage-pipeline-02-view.md`, Task 2).
///
/// The continuum driver syncs `chain.coords` into `sys.xcart` at the END of a
/// round only; the abort completion loop then keeps moving atoms with
/// `force_place` without syncing, so today's writeback block must read
/// `chain.coords` to be correct. Moving `xcart` to be the single home (and
/// the writeback to `RigidView::capture_from_xcart`, which reads `ctx.xcart`)
/// must not move a single bit of the answer. These literals pin that.
///
/// Provenance: captured 2026-09-03 from the build at commit c8fb40e
/// (`cargo test -p molcrafts-molpack`, default features, debug profile) with
/// a scratch harness that printed `PackResult::positions()` through `{:?}`
/// (shortest round-trip form) for exactly this fixture. Deterministic by
/// construction: fixed seed 7, fixed stop_after, no wall-clock, no threads.
/// The run aborts after 2 rounds with 6 forced placements, so the abort
/// completion path IS exercised.
#[test]
fn grow_abort_writeback_golden() {
    // Unwrapped lab-frame positions — growth does not wrap into the box, so
    // coordinates outside `[0, 14]` are expected and part of the golden.
    const GOLDEN: [[F; 3]; 18] = [
        [6.960950727955596, 1.2876400376932167, 3.617693024965275],
        [7.68775644288166, 2.457652934396787, 2.951567455497975],
        [9.159888433744518, 2.4482765704772875, 3.3682681697735872],
        [9.40329260787971, 1.3071109923141822, 4.357912313508713],
        [10.814791573019656, 1.4251069244160106, 4.936399678466913],
        [10.885033188836235, 0.675961130655507, 6.268595870886941],
        [11.727187056008203, 5.611401891774847, 12.800862080686782],
        [11.297571317298337, 4.50730764820194, 11.832716317714103],
        [12.51593652550088, 4.009731814063756, 11.052392981703148],
        [13.58912942202339, 5.100014860120379, 11.031419514295132],
        [14.126112885232246, 5.2607651229180465, 9.607794523909426],
        [15.64827655186169, 5.106309479190209, 9.615619235650529],
        [2.091397543916382, 10.61917964641629, 5.627510535377444],
        [1.4579490039634453, 11.497973776561347, 4.547064412679006],
        [-0.05835140097951674, 11.293687410805898, 4.546661006877859],
        [-0.40194945716705366, 10.017432296902328, 5.317383563975504],
        [-1.9205392935375896, 9.902896032193846, 5.4645785127296564],
        [-2.2935937549473175, 8.457418264951658, 5.799726828003418],
    ];

    let events = Arc::new(Mutex::new(Vec::new()));
    let (n_beads, copies, bond, box_len) = (6usize, 3usize, 1.53 as F, 14.0 as F);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [box_len; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: 2,
        }))
        .run(&[target], 50)
        .expect("a should_stop abort still returns Ok");

    assert!(
        !result.converged,
        "an aborted grow pack must not be reported as converged"
    );
    assert!(
        result.softened > 0,
        "the golden is only meaningful when the abort completion path ran \
         (force_place counts as softening); got softened = {}",
        result.softened
    );
    let rounds = events.lock().expect("recorder mutex").len();
    assert!(
        rounds < 10,
        "growth recorded {rounds} rounds — this golden pins an ABORT, not a \
         completed schedule"
    );

    let pos = result.positions();
    assert_eq!(
        pos.len(),
        GOLDEN.len(),
        "fixture shape changed: {} atoms vs {} golden rows",
        pos.len(),
        GOLDEN.len()
    );
    for (i, (got, want)) in pos.iter().zip(GOLDEN.iter()).enumerate() {
        for k in 0..3 {
            assert_eq!(
                got[k].to_bits(),
                want[k].to_bits(),
                "atom {i} component {k}: abort writeback moved a bit \
                 (got {}, golden {})",
                got[k],
                want[k]
            );
        }
    }
}

// ── Section: push-off — the explicit free-target chain (门槛 2) ────────────
//
// When growth ends unconverged (softened > 0), the entry says so and stops.
// The rigid push-off is the user-explicit chain (placement-seeding spec):
// the SAME free targets go to `GenCanPack::seeded_from(&grown)`, whose
// phases continue on the coor/x growth wrote (Auhl slow push-off /
// Theodorou–Suter staged relaxation, spec §5.4/§5.7). The seeded run must
// (i) NOT run `initial()` — that re-randomizes every COM/Euler and
// teleports the grown chains before descent even starts — and (ii) run the
// phases with movebad disabled, so molecules move by rigid-body descent
// only.

/// Probe for the seeded hand-off: `init_xcart` is the snapshot at
/// `on_initialized` — the state the GENCAN phases start descending from
/// (`on_initialized` fires inside the seeded solver, AFTER any rogue
/// `initial()` would have scrambled the state). `init_calls` counts runs.
struct HandoffProbe {
    init_xcart: Arc<Mutex<Option<Vec<[F; 3]>>>>,
    init_calls: Arc<Mutex<usize>>,
}

impl Handler for HandoffProbe {
    fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {}

    fn on_initialized(&mut self, sys: &PackContext) {
        *self.init_calls.lock().expect("probe mutex") += 1;
        *self.init_xcart.lock().expect("probe mutex") = Some(sys.xcart.clone());
    }
}

/// The `grow_restraint_infeasible_is_not_silent` setup, reused as THE
/// push-off trigger: a radius-1.0 sphere cannot hold a 6-bead chain, so
/// growth force-places and ends with softened > 0 — the state the seeded
/// chain continues from.
fn infeasible_sphere_target() -> Target {
    Target::new(chain_frame(6, 1.53, true), 2)
        .with_restraint(InsideSphereRestraint::new([15.0, 15.0, 15.0], 1.0))
}

/// Grow the infeasible cell honestly (no hidden continuation).
fn grow_infeasible(l: F, max_loops: usize) -> PackResult {
    CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .run(&[infeasible_sphere_target()], max_loops)
        .expect("the infeasible-restraint grow reports through converged/softened, not Err")
}

/// THE push-off no-teleport test (numerics-review regression guard), on the
/// explicit chain.
///
/// Detectors, sharpest first — each names the regression it catches:
///
/// 1. xcart at the seeded run's `on_initialized` is BITWISE the grown
///    result's positions (zero-conversion chaining: the seed carries
///    (coor, rigid) verbatim, and `RigidView::write_xcart` — called from the
///    push-off branch in `gencan/solver.rs` — re-materializes the same
///    xcart from the same bits). A seeded run that re-ran `initial()`
///    would re-randomize every COM/Euler before `on_initialized` fires and
///    break this by Ångströms; a chain that rebuilt the seed from the
///    frame would recompute COMs and break it in the last ulp.
/// 2. Per-copy bonded geometry still equals the template: rigid-body
///    motion is the ONLY move the seeded stage may make on the grown
///    conformers.
/// 3. Honest verdicts on BOTH links: the grow result keeps
///    `converged == false` / `softened > 0`; the seeded GENCAN stage
///    reports its own outcome (`softened == 0` always on the rigid path).
///
/// Deliberately NOT asserted: a bound on how far molecules travel during
/// the push-off. In this setup the travel is legitimately large (~19 Å:
/// force-placed chains get dragged to the sphere by restraint descent), so
/// distance cannot separate descent from a movebad teleport here — and with
/// 2 molecules at `perturb_fraction` 0.05 a rogue movebad would move
/// ⌊0.05·2⌋ = 0 molecules anyway. Detector 1 is the load-bearing one.
#[test]
fn free_chain_push_off_starts_from_grown_state() {
    let (center, radius, l, n_beads, copies) =
        ([15.0 as F, 15.0, 15.0], 1.0 as F, 30.0 as F, 6usize, 2usize);
    let grown = grow_infeasible(l, 3);
    assert!(
        !grown.converged,
        "an unsatisfiable restraint must not report converged"
    );
    assert!(
        grown.softened > 0,
        "softened = 0 — without softening this test has no trigger"
    );

    let init_xcart = Arc::new(Mutex::new(None));
    let init_calls = Arc::new(Mutex::new(0usize));
    let pushed = GenCanPack::new()
        .seeded_from(&grown)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_handler(Box::new(HandoffProbe {
            init_xcart: Arc::clone(&init_xcart),
            init_calls: Arc::clone(&init_calls),
        }))
        .run(&[infeasible_sphere_target()], 3)
        .expect("the seeded push-off returns Ok");

    // Detector 3: honest verdicts on both links.
    assert_eq!(pushed.softened, 0, "GENCAN reports its own softened count");
    let pos = pushed.positions();
    assert_eq!(pos.len(), copies * n_beads, "2 copies × 6 beads");
    let max_center_dist = pos.iter().map(|&p| vdist(p, center)).fold(0.0, F::max);
    assert!(
        max_center_dist > radius + 1e-6,
        "every atom within the radius-1 sphere — the restraint was \
         satisfiable after all and the push-off trigger is gone"
    );

    // Detector 1: exactly one hand-off, starting bitwise from the grown
    // coordinates.
    assert_eq!(
        *init_calls.lock().expect("probe mutex"),
        1,
        "the seeded run initializes exactly once"
    );
    let ginit = init_xcart
        .lock()
        .expect("probe mutex")
        .clone()
        .expect("init_calls == 1 guarantees this");
    let gpos = grown.positions();
    assert_eq!(ginit.len(), gpos.len(), "same xcart length");
    for (i, (a, b)) in ginit.iter().zip(&gpos).enumerate() {
        for k in 0..3 {
            assert_eq!(
                a[k].to_bits(),
                b[k].to_bits(),
                "atom {i} axis {k}: {} vs {} — the seeded run must start \
                 BITWISE from the grown coordinates (zero-conversion \
                 chaining, placement-seeding spec)",
                a[k],
                b[k]
            );
        }
    }

    // Detector 2: the seeded stage may move copies only rigidly.
    let template = zigzag_coords(n_beads, 1.53);
    for c in 0..copies {
        let chain = &pos[c * n_beads..(c + 1) * n_beads];
        for i in 0..n_beads - 1 {
            let want = vdist(template[i], template[i + 1]);
            let d = min_image_cubic(vsub(chain[i + 1], chain[i]), l);
            let got = vdot(d, d).sqrt();
            assert!(
                (got - want).abs() < 1e-6,
                "copy {c} bond ({i},{j}) = {got} vs template {want} — the \
                 grown conformer must survive the push-off rigidly",
                j = i + 1
            );
        }
    }
}

/// Push-off determinism on the chain: grow + seeded GENCAN, twice, must be
/// bitwise reproducible under fixed seeds. A seeded run that consumed an
/// unseeded RNG — e.g. by re-running `initial()` or re-enabling
/// perturbations — or iterated an unordered container would diverge.
#[test]
fn free_chain_push_off_deterministic() {
    let chain = || {
        let grown = grow_infeasible(30.0, 3);
        assert!(grown.softened > 0, "the cell must exercise the push-off");
        GenCanPack::new()
            .seeded_from(&grown)
            .with_seed(7)
            .with_tolerance(2.0)
            .run(&[infeasible_sphere_target()], 3)
            .expect("the seeded push-off returns Ok")
    };
    let a = chain();
    let b = chain();
    assert_eq!(
        a.converged, b.converged,
        "identical chains must reach the same verdict"
    );
    let (pa, pb) = (a.positions(), b.positions());
    assert_eq!(pa.len(), pb.len(), "same atom count");
    for (i, (p, q)) in pa.iter().zip(&pb).enumerate() {
        for k in 0..3 {
            assert_eq!(
                p[k].to_bits(),
                q[k].to_bits(),
                "atom {i} axis {k}: {} vs {} — the seeded chain must be \
                 bitwise reproducible under fixed seeds",
                p[k],
                q[k]
            );
        }
    }
}

// ── Section: Tasks 8-10 — density-resolved box, CG angle prior ─────────────
//
// Task 8 (spec §7, ac-005): `with_density(rho)` on the shared engine
// settings resolves in stage ① to a CUBIC periodic box `[0, L]³` (all axes
// periodic) with `L = cbrt(total_mass_amu / (N_A · rho) · 1e24)` Å (rho in
// g/cm³), the total mass summing over ALL targets × their counts. Masses
// default to element lookup; `Target::with_mass(amu)` overrides the per-copy
// total (the only route for element-"X" targets). Named errors:
// `PackError::DensityConflictsWithBox` (density + explicit box/cell) and
// `PackError::UnknownMass { target }` (density given, a target's mass
// unresolvable, no override). Density is solver-agnostic: it belongs to the
// shared engine settings, not to `GrowConfig`.
//
// Task 9's in-pack Grow+Gencan composition is gone with the monolithic
// entry (engine-entry-split): the explicit chain over a fixed matrix is
// `grow_then_gencan_chaining_over_fixed_matrix` below.
//
// Task 10, CG half (spec §4a′/§5.5, ac-011): `AnglePrior` makes the bond
// angle a sampling degree of freedom. `AnglePrior::Template` (the default)
// copies template angles verbatim — the AA behavior; `AnglePrior::Wlc`
// (via `wlc_from_c_inf`) is the CG path, calibrated so a discrete worm-like
// chain reproduces c∞ = (1+⟨cosθ′⟩)/(1−⟨cosθ′⟩), ⟨cosθ′⟩ = (c∞−1)/(c∞+1),
// where θ′ is the bond-deflection angle.

// The Task 10 CG surface (compile-RED until `AnglePrior` lands in
// `src/grow/prior.rs` and is re-exported from `molpack::grow`).
use molpack::grow::AnglePrior;

/// CODATA Avogadro constant, exactly as fixed by the 2019 SI.
const AVOGADRO: F = 6.02214076e23;

/// The frozen stage-① density formula (spec §7, ac-005):
/// `L = cbrt(total_mass_amu / (N_A · rho) · 1e24)` Å, rho in g/cm³.
fn expected_cubic_edge(total_mass_amu: F, rho: F) -> F {
    (total_mass_amu / (AVOGADRO * rho) * 1e24).cbrt()
}

/// Bond-graph unwrap of a linear `(i, i+1)` chain under the minimum image in
/// a cubic periodic box of edge `l`: the sequential walk that
/// [`linear_chain_rg`] performs, exposed so end-to-end vectors and interior
/// angles can be measured on chains that wrap the box.
fn unwrap_linear_chain(xyz: &[[F; 3]], l: F) -> Vec<[F; 3]> {
    let mut un: Vec<[F; 3]> = Vec::with_capacity(xyz.len());
    un.push(xyz[0]);
    for w in xyz.windows(2) {
        let prev = *un.last().expect("unwrap walk is never empty");
        let d = min_image_cubic(vsub(w[1], w[0]), l);
        un.push([prev[0] + d[0], prev[1] + d[1], prev[2] + d[2]]);
    }
    un
}

/// Task 8 headline (ac-005): `with_density` resolves to the cubic periodic
/// box `[0, L]³` with L from the frozen formula, the output frame carries
/// that box, and the density recomputed from (total mass, box volume) hits
/// the declared rho to 1e-6 relative — the analytic round-trip that catches
/// unit mistakes.
#[test]
fn density_resolves_cubic_box() {
    let rho = 0.05 as F;
    let mass_per_copy = 72.0 as F; // Target::with_mass — the helper frame has
    let copies = 2usize; // no element column, so masses are "X"
    let target = Target::new(chain_frame(5, 1.5, true), copies).with_mass(mass_per_copy);
    let result = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(rho)
        .run(&[target], 50)
        .expect(
            "a Grow pack at rho = 0.05 g/cm³ is roomy and must succeed — the density IS the box",
        );

    let total_mass = mass_per_copy * copies as F;
    let l = expected_cubic_edge(total_mass, rho);
    let simbox =
        result.frame.simbox.as_ref().expect(
            "with_density must stamp the output frame's simbox (ac-005 measures its volume)",
        );
    assert_eq!(
        simbox.pbc(),
        [true; 3],
        "the density-resolved box is periodic on ALL axes (spec §7)"
    );
    let lengths = simbox.lengths();
    let origin = simbox.origin_view();
    let mut volume = 1.0 as F;
    for k in 0..3 {
        assert!(
            ((lengths[k] - l) / l).abs() < 1e-9,
            "axis {k}: box edge = {}, expected L = {l} within 1e-9 relative — \
             L = cbrt(M_amu / (N_A·rho) · 1e24) over ALL targets × counts",
            lengths[k]
        );
        assert!(
            origin[k].abs() < 1e-9,
            "axis {k}: origin = {}, the density-resolved box is [0, L]³",
            origin[k]
        );
        volume *= lengths[k];
    }
    let rho_actual = total_mass / AVOGADRO / (volume * 1e-24);
    assert!(
        ((rho_actual - rho) / rho).abs() < 1e-6,
        "mass density recomputed from the output box = {rho_actual} g/cm³, \
         declared rho = {rho} — must match within 1e-6 relative (ac-005)"
    );
}

/// `with_density` together with an explicit box is a named
/// `DensityConflictsWithBox` error — one source of truth for the volume,
/// never a silent precedence rule (spec §7).
#[test]
fn density_conflicts_with_box() {
    let target = Target::new(chain_frame(5, 1.5, true), 2).with_mass(72.0);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(0.05)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[target], 20)
        .expect_err("with_density combined with with_periodic_box must be rejected");
    assert!(
        matches!(err, PackError::DensityConflictsWithBox),
        "expected PackError::DensityConflictsWithBox, got {err:?}"
    );
}

/// `with_density` on a target whose element masses cannot be resolved
/// (element "X", no `with_mass` override) is a named
/// `UnknownMass { target }` error — molpack does not guess masses (spec §7,
/// principle 3).
#[test]
fn density_needs_mass() {
    // chain_frame has no element column → Target::new resolves "X" masses.
    let target = Target::new(chain_frame(5, 1.5, true), 2);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(0.05)
        .run(&[target], 20)
        .expect_err("density with an unresolvable target mass and no with_mass must be rejected");
    assert!(
        matches!(err, PackError::UnknownMass { target: 0 }),
        "expected PackError::UnknownMass {{ target: 0 }}, got {err:?}"
    );
}

/// Density is solver-agnostic (spec §7): it lives on the shared engine
/// settings, stage ①, and the rigid-body `GenCanPack` entry resolves the
/// same cubic box from the same formula.
#[test]
fn density_works_for_gencan() {
    let rho = 0.02 as F;
    let target = Target::new(chain_frame(5, 1.5, true), 2).with_mass(72.0);
    let result = GenCanPack::new()
        .with_seed(3)
        .with_tolerance(2.0)
        .with_density(rho)
        .run(&[target], 20)
        .expect(
            "a pure-Gencan pack with with_density must succeed — density is not a Grow feature",
        );

    let l = expected_cubic_edge(2.0 * 72.0, rho);
    let simbox = result
        .frame
        .simbox
        .as_ref()
        .expect("the gencan path must stamp the density-resolved simbox too");
    assert_eq!(simbox.pbc(), [true; 3], "cubic periodic box on all axes");
    let lengths = simbox.lengths();
    for k in 0..3 {
        assert!(
            ((lengths[k] - l) / l).abs() < 1e-9,
            "axis {k}: box edge = {}, expected L = {l} within 1e-9 relative \
             — same stage-① formula regardless of solver",
            lengths[k]
        );
    }
}

/// THE ac-011 test: 150 programmatic Kremer–Grest bead chains (100 beads,
/// bond 0.97σ, declared contact 0.85σ via tolerance 0.85, `exclusion_depth`
/// explicitly 2), grown at the KG MELT number density ρ* = 0.85σ⁻³ with
/// uniform torsions and a WLC angle prior calibrated to c∞ = 1.76, measure
/// back c_n = 1.76 ± 10%. The interior-angle spread across copies proves the
/// angle actually became a sampling degree of freedom
/// (`AnglePrior != Template` path exercised).
///
/// Why the MELT ensemble and not isolated chains: with 1-4-and-beyond
/// hard-core self-avoidance, an isolated grown chain is a
/// (Rosenbluth-biased) self-avoiding walk and can never measure 1.76 —
/// instrumented runs gave c_n = 2.80 (depth 2) and 2.50 (depth 3) in a
/// dilute 300 Å box. Auhl's c∞ = 1.76 is the melt value: interchain
/// crowding screens the intrachain excluded volume (spec §5.5's
/// "Kremer–Grest melt" anchor). The spec's ladder step 2 (isolated chain +
/// hard core: measure and REPORT, no assertion) and ac-011 (assert 1.76)
/// are consistent exactly because the ac-011 measurement is a melt one.
/// Asserting here bets that growth AT density reproduces the screening —
/// the same scientific bet ac-006(3) makes for the AA melt.
///
/// Why tolerance 0.85σ and not the full bead diameter 1.0σ: a strict 1.0σ
/// hard core at ρ* = 0.85 is hard-sphere packing fraction η ≈ 0.445 — KG
/// melts only exist there because WCA is soft, and an instrumented run at
/// tolerance 1.0 degraded accordingly (softened = 1372, final fdist ≈ 1.0:
/// forced-placement leftovers). The packer's tolerance is a
/// PRE-RELAXATION contact criterion, not the interaction diameter —
/// Packmol's canonical 2.0 Å is ~0.6-0.7 of the heavy-atom σ. The KG
/// equivalent is 0.85σ: it matches the Auhl slow-push-off floor of 0.8σ
/// that `min_hard_scale`'s default cites (spec §4f / §5.4, the same
/// physics) and brings η to 0.445 × 0.85³ ≈ 0.27, comfortably feasible.
/// The constructive guarantee stays strict — min separation ≥ the declared
/// tolerance, bitwise fdist == 0, softened == 0 — it is the declared
/// contact that changes to the physically meaningful one.
///
/// Why depth 2 and not 1: ac-011's claim is that the angle-prior
/// CALIBRATION hits c∞, so the 1-3 pair must be excluded like the 1-2 —
/// the 1-3 distance is a pure function of the bond angle
/// (d13 = 2·b·sin(θ/2)), so hard-core scoring it truncates the WLC
/// deflection distribution (at the then-declared contact 1.0σ vs bond
/// 0.97σ it blocked θ' ≳ 118°) and stacks hard-core stiffening ON TOP of
/// the prior: measured c_n = 3.02 at depth 1 (and 4.70 before the field
/// stopped soft-shell-charging same-molecule pairs). At depth 2 the
/// 1-4-and-beyond pairs stay hard-core scored — the chain still cannot
/// thread itself. A shallower depth is a different physical model and
/// would need an empirical recalibration of κ, not a wider band.
#[test]
fn grow_cg_kremer_grest_c_inf() {
    let (n_beads, copies) = (100usize, 150usize);
    let bond = 0.97 as F;
    let tolerance = 0.85 as F;
    // Default `min_hard_scale` (GrowConfig::new) — the Auhl push-off floor.
    let min_hard_scale = 0.8 as F;
    // KG melt number density ρ* = 0.85σ⁻³, σ = 1.0 (the declared contact is
    // 0.85σ — see the doc comment): L = (N / ρ*)^(1/3) ≈ 26.03σ for
    // N = 15000 beads.
    let l = ((n_beads * copies) as F / 0.85).cbrt();
    let events: StepEvents = Arc::new(Mutex::new(Vec::new()));
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_angle_prior(AnglePrior::wlc_from_c_inf(1.76))
        .with_exclusion_depth(2);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);
    let result = CbmcGrow::from_config(cfg)
        .with_seed(5)
        .with_tolerance(tolerance)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: usize::MAX,
        }))
        .run(&[target], 200)
        .expect("the KG melt at ρ* = 0.85σ⁻³ must grow — this is the CG headline pack (ac-011)");

    let pos = result.positions();
    assert_eq!(pos.len(), copies * n_beads, "150 copies × 100 beads");
    let mol_of: Vec<usize> = (0..pos.len()).map(|i| i / n_beads).collect();

    // ── What the packed state must satisfy ───────────────────────────────
    //
    // TODO(grow-axes ac-004): restore the strict fdist == 0 assertion once
    // softening is per-chain and recoverable.
    let floor = min_hard_scale * tolerance;
    if result.softened == 0 {
        assert_eq!(
            result.fdist.to_bits(),
            (0.0 as F).to_bits(),
            "fdist = {} with softened == 0 — an UNSOFTENED growth run keeps \
             the strict constructive hard-core guarantee (ac-004): candidates \
             are rejected, never penalized",
            result.fdist
        );
    } else {
        // Softened: the guarantee weakens from `tolerance` to the softening
        // floor, but it stays CONSTRUCTIVE — every committed candidate was
        // still hard-rejected against `hard_scale × tolerance ≥ floor`, so no
        // pair may sit below the floor. Two independent rulers: the O(N²)
        // worst pair and the cell-binned population count.
        let dmin = min_inter_distance(&pos, &mol_of, l);
        assert!(
            dmin >= floor - 1e-9,
            "softened = {} and the closest inter-molecular pair is {dmin} — \
             below the softening floor {floor} (= min_hard_scale {min_hard_scale} \
             × tolerance {tolerance}). Softening lowers the declared contact; \
             it must never abandon hard rejection",
            result.softened
        );
        assert_eq!(
            count_inter_pairs_below(&pos, &mol_of, l, floor - 1e-9),
            0,
            "inter-molecular pairs below the softening floor {floor} — the \
             floor is a bound on the whole population, not just on the \
             worst pair (min = {dmin})"
        );
    }

    // ── How softening is allowed to happen: a ladder, never a free-fall ──
    let events = events.lock().expect("recorder mutex");
    assert_radscale_ladder(&events, min_hard_scale, "KG melt at ρ* = 0.85");

    // ── The science: melt statistics of the grown chains ─────────────────
    let mut sum_r2 = 0.0 as F;
    let mut angle_lo = F::INFINITY;
    let mut angle_hi = F::NEG_INFINITY;
    for c in 0..copies {
        let un = unwrap_linear_chain(&pos[c * n_beads..(c + 1) * n_beads], l);
        let r = vsub(un[n_beads - 1], un[0]);
        sum_r2 += vdot(r, r);
        let a = vangle(un[49], un[50], un[51]);
        angle_lo = angle_lo.min(a);
        angle_hi = angle_hi.max(a);
    }
    let c_n = sum_r2 / copies as F / ((n_beads - 1) as F * bond * bond);
    assert!(
        (c_n - 1.76).abs() <= 0.176,
        "KG chains: c_n = {c_n:.4}, target c∞ = 1.76 ± 10% (Auhl et al. \
         NRRW value, ac-011)"
    );
    assert!(
        angle_hi - angle_lo > 0.1,
        "interior angle at bead 50 spans [{angle_lo:.4}, {angle_hi:.4}] rad \
         across copies — with AnglePrior::Template every copy would show the \
         template angle; the Wlc path must actually sample angles (ac-011)"
    );
}

/// The round loop must be BOUNDED — a configuration the strict hard core
/// cannot satisfy has to END, not spin.
///
/// Fixture: 20 × 24-bead chains (bond 1.53) in a 17 Å periodic box at
/// tolerance 2.0 — reduced density N·tol³/L³ ≈ 0.78 — with the softening
/// escape switched OFF from both sides: `min_hard_scale(1.0)` forbids any
/// shrink and `soften_after(usize::MAX / 4)` puts the per-chain rung (and
/// with it the `2 × soften_after` forced-placement escape) out of reach. The
/// growth is then geometrically impossible AND has no valve, which is exactly
/// the state the round loop must survive.
///
/// Why 17 Å: the same fixture converges constructively at 20 Å (0.24 s) and
/// at 19 Å (8.8 s), and livelocks at ≤ 18 Å. 17 Å is inside the livelock
/// regime with margin (measured: 32 974 rounds in 25 s and still running),
/// while staying small enough that the fix's cap — `max_loops.max(1) ×
/// (n_steps + 1)` = 3 × 22 = 66 rounds at ~0.5 ms per round — is reached in
/// milliseconds.
///
/// Contract on the way out: `Ok`, `converged == false` (nothing was proved),
/// `softened > 0` (the forced completions are counted like every other
/// break of the constructive guarantee), finite coordinates, and — as in
/// [`grow_abort_keeps_bonded_geometry`] — bonded geometry still chemical,
/// because forced placement builds from the template's internal coordinates
/// (measured worst deviation on a completing run: 6e-15 Å).
#[test]
fn grow_dense_strict_core_terminates_unconverged() {
    let (copies, n_beads, bond, l) = (20usize, 24usize, 1.53 as F, 17.0 as F);
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_min_hard_scale(1.0)
        .with_soften_after(usize::MAX / 4);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);

    // The run is the thing under test, so it may not be allowed to hang the
    // suite: it gets its own thread and a wall-clock budget two orders of
    // magnitude above the capped cost.
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let outcome = CbmcGrow::from_config(cfg)
            .with_seed(11)
            .with_tolerance(2.0)
            .with_periodic_box([0.0; 3], [l; 3], [true; 3])
            .run(&[target], 3);
        let _ = tx.send(outcome);
    });
    let result = match rx.recv_timeout(std::time::Duration::from_secs(120)) {
        Ok(outcome) => outcome.expect(
            "a growth run that cannot satisfy its hard core must still return \
             Ok with converged == false — an impossible density is a result, \
             not an error",
        ),
        Err(std::sync::mpsc::RecvTimeoutError::Timeout) => panic!(
            "growth did not return within 120 s: 20 × 24 beads in a 17 Å box \
             at tolerance 2.0 with a strict hard core (min_hard_scale = 1.0) \
             cannot be satisfied and the round loop (src/grow/driver.rs) has \
             no upper bound, so it spins forever. The driver must cap the \
             round loop and force-complete the pending chains through the \
             abort path (debt D-01 (ii))"
        ),
        Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
            panic!("the growth thread died without sending a result")
        }
    };

    assert!(
        !result.converged,
        "a run that had to force its way out of an unsatisfiable hard core \
         must report converged == false — the cap is a surrender, not a proof"
    );
    assert!(
        result.softened > 0,
        "softened = 0 after a capped run — every forced completion breaks the \
         constructive guarantee and must be counted, or the caller cannot \
         tell a capped result from a clean one"
    );

    let pos = result.positions();
    assert_eq!(pos.len(), copies * n_beads, "20 copies × 24 beads");
    for (i, p) in pos.iter().enumerate() {
        assert!(
            p.iter().all(|v| v.is_finite()),
            "atom {i} at {p:?} — a capped run must still write real \
             coordinates for every atom"
        );
    }

    let bonds = result
        .frame
        .get("bonds")
        .expect("assembled frame keeps bonds");
    let ai = bonds.get_uint("atomi").expect("atomi");
    let aj = bonds.get_uint("atomj").expect("atomj");
    assert_eq!(
        ai.len(),
        copies * (n_beads - 1),
        "template 1-2 bonds must be tiled onto every copy"
    );
    let mut lo = F::INFINITY;
    let mut hi = F::NEG_INFINITY;
    for (&a, &b) in ai.iter().zip(aj.iter()) {
        let d = min_image_cubic(vsub(pos[a as usize], pos[b as usize]), l);
        let len = vdot(d, d).sqrt();
        lo = lo.min(len);
        hi = hi.max(len);
    }
    assert!(
        (lo - bond).abs() <= 1e-6 && (hi - bond).abs() <= 1e-6,
        "bonded range [{lo}, {hi}] Å must stay the template bond {bond} Å to \
         1e-6 — forced completion places atoms from the template's internal \
         coordinates, so it cannot invent 0-length or box-scale bonds"
    );
}

/// Softening must be EARNED by repeated dead ends — running out of budget is
/// not a licence to shrink the hard core.
///
/// Historical framing: before the D-01 fix the driver carried a
/// `regrow_budget` (`max_loops × n_chains` regrow events), and once that
/// budget was exhausted every further dead end shrank the hard core directly.
/// That path is gone — the driver now terminates through the round cap — and
/// `max_loops = 1` is kept here as the cheap witness that the replacement still
/// routes every softening rung through the per-chain `soften_after` ladder.
///
/// Same KG melt fixture as [`grow_cg_kremer_grest_c_inf`], run with
/// `max_loops = 1` — which is what made the old budget run out long before 150
/// chains of 100 beads were grown.
///
/// The observables, all consequences of one commit attempt per chain per
/// round: growth needs at least `n_beads - 3` rounds; no rung may fall before
/// round `soften_after`, because that many dead ends cannot have accumulated
/// on any chain yet; and every rung is a single 3 % step. Before the fix the
/// first rung fell at round 69 and took the scale 1.0 → 0.8 inside that single
/// round — exactly what the ladder assertion below rejects.
#[test]
fn grow_softening_needs_repeated_dead_ends() {
    let (n_beads, copies) = (100usize, 150usize);
    let bond = 0.97 as F;
    // GrowConfig::new defaults, spelled out because the assertions below are
    // arithmetic on them.
    let (soften_after, min_hard_scale) = (50usize, 0.8 as F);
    let l = ((n_beads * copies) as F / 0.85).cbrt();
    let events: StepEvents = Arc::new(Mutex::new(Vec::new()));
    let cfg = GrowConfig::new(TorsionPrior::Uniform)
        .with_angle_prior(AnglePrior::wlc_from_c_inf(1.76))
        .with_exclusion_depth(2);
    let target = Target::new(chain_frame(n_beads, bond, true), copies);
    let result = CbmcGrow::from_config(cfg)
        .with_seed(5)
        .with_tolerance(0.85)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .with_handler(Box::new(Recorder {
            events: Arc::clone(&events),
            stop_after: usize::MAX,
        }))
        .run(&[target], 1)
        .expect("a max_loops = 1 run is not an error — the run must return Ok");
    assert_eq!(
        result.positions().len(),
        copies * n_beads,
        "150 copies × 100 beads"
    );

    let events = events.lock().expect("recorder mutex");
    // A 100-bead chain is a 3-atom seed plus 97 steps, one commit per round.
    assert!(
        events.len() >= n_beads - 3,
        "growth recorded {} rounds for a {n_beads}-bead chain — it cannot \
         have grown anything in fewer than {} rounds, so this fixture is not \
         exercising the budget path at all",
        events.len(),
        n_beads - 3
    );

    let first_rung = events.iter().find(|&&(_, radscale, _, _)| radscale < 1.0);
    if let Some(&(round, radscale, _, _)) = first_rung {
        assert!(
            round >= soften_after,
            "the hard core first softened to {radscale} at round {round}, \
             before any chain could have reached {soften_after} dead ends (a \
             chain dead-ends at most once per round). Every rung must be \
             earned on the per-chain soften_after ladder; the now-removed \
             regrow_budget path used to shrink the core on a chain's first \
             dead end once its budget ran out (debt D-01 (i))"
        );
    }
    assert_radscale_ladder(&events, min_hard_scale, "KG melt, max_loops = 1");
}

/// `AnglePrior::Template` is the default AND the AA behavior: a default
/// `GrowConfig` and one with the explicit `Template` prior produce
/// bitwise-identical packs, and the constructive `fdist == 0.0` guarantee
/// (ac-004) is untouched by the angle-prior seam.
#[test]
fn angle_prior_template_is_default() {
    let dbg = format!("{:?}", GrowConfig::new(TorsionPrior::Uniform));
    assert!(!dbg.is_empty(), "GrowConfig must keep a useful Debug impl");

    let pack = |cfg: GrowConfig| {
        CbmcGrow::from_config(cfg)
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
            .run(&[Target::new(chain_frame(8, 1.53, true), 4)], 50)
            .expect("a small grow pack in a roomy box must succeed")
    };

    let default_result = pack(GrowConfig::new(TorsionPrior::Uniform));
    let explicit_result =
        pack(GrowConfig::new(TorsionPrior::Uniform).with_angle_prior(AnglePrior::Template));

    for result in [&default_result, &explicit_result] {
        assert_eq!(
            result.fdist.to_bits(),
            (0.0 as F).to_bits(),
            "fdist = {} — the ac-004 constructive guarantee must survive the \
             angle-prior seam",
            result.fdist
        );
        assert_eq!(result.softened, 0, "roomy box: no softening");
    }

    let (pa, pb) = (default_result.positions(), explicit_result.positions());
    assert_eq!(pa.len(), pb.len(), "same atom count");
    for (i, (a, b)) in pa.iter().zip(&pb).enumerate() {
        for k in 0..3 {
            assert_eq!(
                a[k].to_bits(),
                b[k].to_bits(),
                "atom {i} axis {k}: {} vs {} — an explicit AnglePrior::Template \
                 must be bitwise identical to the default, pinning Template as \
                 the AA default (spec §4a′)",
                a[k],
                b[k]
            );
        }
    }
}

/// Cell-occupancy bookkeeping behind void-biased seeding: the empty-cell
/// count tracks insert/remove exactly, and `empty_cell_point` never offers
/// a point inside an occupied cell.
#[test]
fn field_empty_cell_bookkeeping() {
    // 4³ box, cutoff 2 → 2×2×2 grid, cell edge 2.
    let mut field = OverlapField::new(
        [0.0; 3],
        [4.0; 3],
        [true; 3],
        vec![0.5; 4],
        vec![0, 1, 2, 3],
        vec![0; 4],
        2.0,
    );
    assert_eq!(field.n_empty_cells(), 8);

    field.insert(0, [1.0, 1.0, 1.0]); // cell (0,0,0)
    assert_eq!(field.n_empty_cells(), 7);
    field.insert(1, [1.5, 1.5, 1.5]); // same cell
    assert_eq!(field.n_empty_cells(), 7);

    // Every pick must avoid the occupied cell [0,2)³.
    for i in 0..64 {
        let u0 = (i as F + 0.5) / 64.0;
        let p = field
            .empty_cell_point([u0, 0.5, 0.5, 0.5])
            .expect("empty cells exist");
        assert!(
            p.iter().any(|&x| x >= 2.0),
            "sampled point {p:?} landed in the occupied cell"
        );
    }

    field.remove(1);
    assert_eq!(field.n_empty_cells(), 7, "cell still holds slot 0");
    field.remove(0);
    assert_eq!(field.n_empty_cells(), 8);

    // A fully occupied grid offers no cavity.
    let mut full = OverlapField::new(
        [0.0; 3],
        [2.0; 3],
        [true; 3],
        vec![0.5],
        vec![0],
        vec![0],
        2.0,
    );
    full.insert(0, [1.0, 1.0, 1.0]);
    assert_eq!(full.n_empty_cells(), 0);
    assert!(full.empty_cell_point([0.3, 0.5, 0.5, 0.5]).is_none());
}

// ── Section: engine-entry-split — the growth entry ─────────────────────────

/// A ring template is refused by name, never silently grown open (the old
/// tree decomposition dropped ring-closing bonds without a word).
#[test]
fn ring_template_is_refused() {
    let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.75, 1.3, 0.0]];
    let frame = frame_from_parts(&coords, &[(0, 1), (1, 2), (2, 0)]);
    let err = CbmcGrow::new(TorsionPrior::Uniform)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(frame, 2)], 60)
        .expect_err("a ring must be refused");
    let msg = format!("{err}");
    assert!(msg.contains("ring"), "named rejection, got: {msg}");
}

/// The seeded chain's named rejections and composability.
///
/// (Bitwise parity of the free chain against the deleted
/// `CbmcGrow::with_push_off` path was proven by a one-time migration test
/// before that knob was removed — placement-seeding spec, numerical
/// contract.)
#[test]
fn seeded_run_contract() {
    let grown = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(9)
        .with_tolerance(1.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 60)
        .expect("grow stage runs");
    assert!(grown.converged);

    // Shape mismatch is a named rejection, not a scrambled pack.
    let err = GenCanPack::new()
        .seeded_from(&grown)
        .run(&[Target::new(chain_frame(5, 1.5, true), 3)], 10)
        .expect_err("a seed for 2 copies must refuse 3");
    assert!(
        matches!(err, PackError::SeedMismatch { .. }),
        "expected SeedMismatch, got {err:?}"
    );

    // The cell travels with the seed; a second box is the existing
    // mutual-exclusion error, never a silent precedence rule.
    let err = GenCanPack::new()
        .seeded_from(&grown)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 10)
        .expect_err("seed cell + declared box must be rejected");
    let msg = format!("{err}");
    assert!(msg.contains("mutually exclusive"), "got: {msg}");

    // A fixed matrix may be appended after the seeded free targets.
    let dimer = Target::from_coords(&[[0.0; 3], [1.5, 0.0, 0.0]], &[0.5, 0.5], 1)
        .with_centering(molpack::CenteringMode::Off)
        .fixed_at([2.0, 2.0, 2.0]);
    let packed = GenCanPack::new()
        .seeded_from(&grown)
        .with_seed(3)
        .with_tolerance(1.0)
        .run(&[Target::new(chain_frame(5, 1.5, true), 2), dimer], 40)
        .expect("seeded run with an appended fixed matrix runs");
    assert_eq!(packed.natoms(), grown.natoms() + 2);
}

/// Explicit fixed-matrix chaining: grow chains, freeze them via
/// `Target::fixed_from`, and pack rigid molecules around them with the
/// GENCAN entry. The matrix must come back verbatim.
#[test]
fn grow_then_gencan_chaining_over_fixed_matrix() {
    let grown = CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(9)
        .with_tolerance(1.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(chain_frame(5, 1.5, true), 2)], 60)
        .expect("grow stage runs");
    assert!(grown.converged);

    let matrix = Target::fixed_from(&grown);
    let dimer = Target::from_coords(&[[0.0; 3], [1.5, 0.0, 0.0]], &[0.5, 0.5], 4);
    let packed = GenCanPack::new()
        .with_seed(3)
        .with_tolerance(1.0)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[matrix, dimer], 80)
        .expect("chained rigid stage runs");

    assert!(
        packed.converged,
        "rigid stage around the frozen chains converges"
    );
    assert_eq!(packed.natoms(), grown.natoms() + 8);
    // The frozen matrix atoms are the first target and must be unchanged.
    let (a, b) = (grown.positions(), packed.positions());
    for (pa, pb) in a.iter().zip(b.iter()) {
        for k in 0..3 {
            assert_eq!(pa[k].to_bits(), pb[k].to_bits(), "matrix must not move");
        }
    }
}

// ── Section: lattice growth (lattice-growth-phase spec) ────────────────────

/// A CG bead chain decorates essentially onto the lattice (uniform template
/// bond = lattice bond), so the occupancy guard's 2nd-neighbour distance
/// becomes a real constructive guarantee: fdist == 0.0 strict at easy fill.
#[test]
fn lattice_grow_bead_chain_constructive() {
    use molpack::LatticeGrow;
    let (copies, n_beads, l) = (8usize, 12usize, 26.0 as F);
    let result = LatticeGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .run(&[Target::new(chain_frame(n_beads, 1.53, true), copies)], 60)
        .expect("a lattice grow at easy fill runs");
    assert!(
        result.converged,
        "easy fill must converge (fdist = {}, softened = {})",
        result.fdist, result.softened
    );
    assert_eq!(result.softened, 0, "no escapes at easy fill");
    assert_eq!(
        result.fdist.to_bits(),
        (0.0 as F).to_bits(),
        "fdist = {} — the occupancy guard is constructive on a bead chain",
        result.fdist
    );

    let pos = result.positions();
    assert_eq!(pos.len(), copies * n_beads);
    // Decoration keeps the template's bonded geometry exactly.
    for c in 0..copies {
        let chain = &pos[c * n_beads..(c + 1) * n_beads];
        for i in 0..n_beads - 1 {
            let d = min_image_cubic(vsub(chain[i + 1], chain[i]), l);
            let got = vdot(d, d).sqrt();
            assert!(
                (got - 1.53).abs() < 1e-6,
                "copy {c} bond {i}: {got} vs template 1.53"
            );
        }
    }

    // Same seed, bitwise reproducible.
    let again = LatticeGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .run(&[Target::new(chain_frame(n_beads, 1.53, true), copies)], 60)
        .expect("second run");
    for (a, b) in result.positions().iter().zip(again.positions().iter()) {
        for k in 0..3 {
            assert_eq!(a[k].to_bits(), b[k].to_bits(), "bitwise determinism");
        }
    }
}

/// Branched heavy-atom templates are staged: named rejection, no silent
/// degradation (lattice-growth-phase spec, 拓扑范围).
#[test]
fn lattice_grow_rejects_branched() {
    use molpack::LatticeGrow;
    // A 5-atom star: center bonded to 3 arms + one arm extended.
    let coords = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [-1.5, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [3.0, 0.0, 0.0],
    ];
    let frame = frame_from_parts(&coords, &[(0, 1), (0, 2), (0, 3), (1, 4)]);
    let err = LatticeGrow::new(TorsionPrior::Uniform)
        .with_periodic_box([0.0; 3], BOX_MAX, [true; 3])
        .run(&[Target::new(frame, 2)], 60)
        .expect_err("a branched backbone must be refused by name");
    let msg = format!("{err}");
    assert!(msg.contains("branched"), "named rejection, got: {msg}");
}

/// The melt-density pipeline: lattice growth at an occupancy where the
/// continuum solver grinds, then the explicit seeded push-off on whatever
/// contacts decoration left. Honest verdicts on both links.
#[test]
fn lattice_grow_then_seeded_push_off_dense() {
    use molpack::{GenCanPack, LatticeGrow};
    // 20 × 24-bead chains in a 22 Å box: site occupancy ≈ 0.5 — far past
    // the continuum grow ceiling.
    let (copies, n_beads, l) = (20usize, 24usize, 22.0 as F);
    let target = || Target::new(chain_frame(n_beads, 1.53, true), copies);
    let grown = LatticeGrow::new(TorsionPrior::Uniform)
        .with_seed(11)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [l; 3], [true; 3])
        .run(&[target()], 60)
        .expect("dense lattice grow runs");
    assert_eq!(grown.natoms(), copies * n_beads);

    if grown.converged {
        return; // constructive already — nothing to push off
    }
    let pushed = GenCanPack::new()
        .seeded_from(&grown)
        .with_seed(11)
        .with_tolerance(2.0)
        .run(&[target()], 120)
        .expect("seeded push-off runs");
    assert_eq!(pushed.natoms(), copies * n_beads);
    assert!(
        pushed.fdist <= grown.fdist,
        "push-off must not worsen the contacts ({} -> {})",
        grown.fdist,
        pushed.fdist
    );
}

// Bitwise parity of `CbmcGrow` against the deleted monolithic `Molpack`
// path — pure grow AND the explicit push-off chain — was proven
// test-for-test before that surface was removed (engine-entry-split
// migration record). The determinism gates above
// (`grow_deterministic_same_seed`, `grow_push_off_deterministic`) now run
// on the entry directly.

// ── Section: the seam markers the two growth stages declare ───────────────
//
// Owner-side half of acceptance ac-008 (`stage-pipeline-04-stage`): what
// `GrowStage` and `LatticeStage` declare on the stage seam belongs to their
// owner — `tests/stage.rs` verifies the seam's own behaviour with fakes and
// boots no real algorithm. Both stages build their placements from nothing
// and return with every free molecule placed, so both read
// `Placed::None -> Placed::All`.

/// The continuum growth stage requires nothing placed and guarantees
/// everything placed.
#[test]
fn grow_stage_requires_none_guarantees_all() {
    use molpack::grow::driver::GrowStage;
    use molpack::{Placed, Stage};

    let targets = [Target::new(chain_frame(6, 1.53, true), 2)];
    let stage = GrowStage::from_targets(&targets, &GrowConfig::new(TorsionPrior::Uniform), 7)
        .expect("a bead-chain target with a template builds a growth stage");

    assert_eq!(
        stage.requires().placed,
        Placed::None,
        "growth constructs its own placements, so it requires nothing placed"
    );
    assert_eq!(
        stage.guarantees().placed,
        Placed::All,
        "growth returns with every free molecule placed"
    );
}

/// The diamond-lattice growth stage declares the same two markers.
#[test]
fn lattice_stage_requires_none_guarantees_all() {
    use molpack::grow::lattice::LatticeStage;
    use molpack::{LatticeConfig, Placed, Stage};

    let targets = [Target::new(chain_frame(6, 1.53, true), 2)];
    let stage = LatticeStage::from_targets(&targets, &LatticeConfig::new(TorsionPrior::Uniform), 7)
        .expect("a bead-chain target with a template builds a lattice stage");

    assert_eq!(
        stage.requires().placed,
        Placed::None,
        "lattice growth decorates its own placements onto the lattice"
    );
    assert_eq!(
        stage.guarantees().placed,
        Placed::All,
        "lattice growth returns with every free molecule placed"
    );
}

// ── Section: the stage identity carried on StepInfo ───────────────────────

/// One recorded `StageInfo`, copied out of a `StepInfo` inside `on_step`.
type StageEvents = Arc<Mutex<Vec<(usize, usize, &'static str)>>>;

/// Records the stage identity of every `on_step` event.
struct StageRecorder {
    events: StageEvents,
}

impl Handler for StageRecorder {
    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        self.events.lock().expect("stage recorder mutex").push((
            info.stage.index,
            info.stage.total,
            info.stage.name,
        ));
    }
}

/// Every `StepInfo` names the stage that emitted it. A single-stage run is
/// stage 0 of 1, and its name is the stage's own `name()` — `"growth"` for
/// the continuum growth driver.
#[test]
fn step_info_names_the_stage_that_emitted_it() {
    let events: StageEvents = Arc::new(Mutex::new(Vec::new()));
    let target = Target::new(chain_frame(12, 1.53, true), 8);
    CbmcGrow::new(TorsionPrior::Uniform)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0; 3], [26.0; 3], [true; 3])
        .with_handler(Box::new(StageRecorder {
            events: Arc::clone(&events),
        }))
        .run(&[target], 50)
        .expect("the Task 5 small-melt baseline must still succeed");

    let events = events.lock().expect("stage recorder mutex");
    assert!(
        !events.is_empty(),
        "growth emitted no on_step events, so the stage identity on StepInfo \
         is unobservable"
    );
    for &(index, total, name) in events.iter() {
        assert_eq!(index, 0, "a single-stage run reports index 0");
        assert_eq!(total, 1, "a single-stage run reports total 1");
        assert_eq!(
            name, "growth",
            "StepInfo.stage.name must be the emitting stage's own name()"
        );
    }
}
