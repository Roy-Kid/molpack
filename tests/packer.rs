//! End-to-end packer tests: small molecule packing, PBC, multi-target,
//! fixed target, error cases.

use molpack::{
    F, InsideBoxRestraint, InsideSphereRestraint, Molpack, NullHandler, OutsideSphereRestraint,
    PackError, Target,
};

// ── helpers ────────────────────────────────────────────────────────────────

fn water_positions() -> Vec<[F; 3]> {
    vec![
        [0.0, 0.0, 0.0],    // O
        [0.96, 0.0, 0.0],   // H
        [-0.24, 0.93, 0.0], // H
    ]
}

fn water_radii() -> Vec<F> {
    vec![1.52, 1.20, 1.20]
}

fn single_atom() -> (Vec<[F; 3]>, Vec<F>) {
    (vec![[0.0, 0.0, 0.0]], vec![1.0])
}

// ── basic packing ──────────────────────────────────────────────────────────

#[test]
fn pack_single_water_in_box() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 1).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let result = Molpack::new().with_seed(42).pack_with_report(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 3);
}

#[test]
fn pack_returns_frame() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 1).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let frame = Molpack::new()
        .with_seed(42)
        .pack(&[target], 5)
        .expect("pack should return a frame");
    assert_eq!(frame.get("atoms").and_then(|b| b.nrows()), Some(3));
}

#[test]
fn pack_three_waters_in_box() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 3).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let result = Molpack::new().with_seed(42).pack_with_report(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 9);
}

#[test]
fn pack_in_sphere() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 3)
        .with_restraint(InsideSphereRestraint::new([0.0, 0.0, 0.0], 20.0));
    let result = Molpack::new().with_seed(42).pack_with_report(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 9);
}

// ── multi-target ───────────────────────────────────────────────────────────

#[test]
fn pack_two_targets_same_box() {
    let (coords, radii) = single_atom();
    let t1 = Target::from_coords(&coords, &radii, 3)
        .with_name("A")
        .with_restraint(InsideBoxRestraint::new(
            [0.0, 0.0, 0.0],
            [40.0, 40.0, 40.0],
            [false; 3],
        ));
    let t2 = Target::from_coords(&coords, &radii, 2)
        .with_name("B")
        .with_restraint(InsideBoxRestraint::new(
            [0.0, 0.0, 0.0],
            [40.0, 40.0, 40.0],
            [false; 3],
        ));
    let result = Molpack::new().with_seed(42).pack_with_report(&[t1, t2], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 5);
}

#[test]
fn pack_mixed_free_and_fixed() {
    let (coords, radii) = single_atom();
    let free = Target::from_coords(&coords, &radii, 2).with_restraint(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [20.0, 20.0, 20.0],
        [false; 3],
    ));
    let fixed = Target::from_coords(&coords, &radii, 1).fixed_at([10.0, 10.0, 10.0]);
    let result = Molpack::new()
        .with_seed(42)
        .pack_with_report(&[free, fixed], 5)
        .expect("should succeed");
    assert_eq!(result.natoms(), 3);
    let fixed_pos = result.positions()[2];
    assert!((fixed_pos[0] - 10.0).abs() < 1e-6);
    assert!((fixed_pos[1] - 10.0).abs() < 1e-6);
}

// ── PBC ────────────────────────────────────────────────────────────────────

#[test]
fn pbc_box_packing() {
    // A fully periodic InsideBoxRestraint both confines atoms softly and
    // declares the system PBC; the packer derives PBC from the restraint
    // — no separate with_periodic on Molpack.
    let target = Target::from_coords(&water_positions(), &water_radii(), 2).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [30.0, 30.0, 30.0], [true; 3]),
    );
    let result = Molpack::new().with_seed(42).pack_with_report(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 6);
}

#[test]
fn invalid_pbc_box_rejected() {
    // Zero-extent y-axis on a periodic restraint must be rejected.
    let target = Target::from_coords(&water_positions(), &water_radii(), 1).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [10.0, 0.0, 10.0], [true; 3]),
    );
    let result = Molpack::new().with_seed(7).pack_with_report(&[target], 5);
    assert!(
        matches!(result, Err(PackError::InvalidPBCBox { .. })),
        "expected InvalidPBCBox, got: {result:?}"
    );
}

#[test]
fn conflicting_periodic_boxes_rejected() {
    // Two targets, each declaring a periodic InsideBoxRestraint with
    // different bounds — must fail at pack() time.
    let t1 = Target::from_coords(&water_positions(), &water_radii(), 1)
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [30.0; 3], [true; 3]));
    let t2 = Target::from_coords(&water_positions(), &water_radii(), 1)
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [40.0; 3], [true; 3]));
    let result = Molpack::new().with_seed(42).pack_with_report(&[t1, t2], 5);
    assert!(
        matches!(result, Err(PackError::ConflictingPeriodicBoxes { .. })),
        "expected ConflictingPeriodicBoxes, got: {result:?}"
    );
}

// ── PBC via packer-level with_periodic_box ─────────────────────────────────

#[test]
fn with_periodic_box_caps_cell_grid() {
    // Regression: without a PBC declaration *and* without an inside
    // restraint, the packer used to fall back to inferring a box from
    // the ±sidemax initial placement, driving the cell grid to ~10⁸
    // cells (42 GB RAM, effectively hanging). `with_periodic_box`
    // bypasses that fallback even when no restraint declares a box.
    let (coords, radii) = single_atom();
    let target = Target::from_coords(&coords, &radii, 4); // no restraints at all
    let result = Molpack::new()
        .with_seed(42)
        .with_periodic_box([0.0; 3], [30.0; 3])
        .pack_with_report(&[target], 2);
    assert!(
        result.is_ok(),
        "expected pack to return quickly under a 30 Å periodic box, got: {result:?}"
    );
    assert_eq!(result.unwrap().natoms(), 4);
}

#[test]
fn with_periodic_box_and_restraint_agree() {
    // Packer-level PBC plus an InsideBoxRestraint declaring exactly
    // the same box is fine — derive_periodic_box sees matching
    // declarations and accepts the run.
    let target = Target::from_coords(&water_positions(), &water_radii(), 2).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [30.0, 30.0, 30.0], [true; 3]),
    );
    let result = Molpack::new()
        .with_seed(42)
        .with_periodic_box([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
        .pack_with_report(&[target], 2);
    assert!(result.is_ok(), "expected agreement, got: {result:?}");
    assert_eq!(result.unwrap().natoms(), 6);
}

#[test]
fn with_periodic_box_conflicts_with_restraint() {
    // Packer-level PBC that disagrees with a restraint-level periodic
    // box on bounds → ConflictingPeriodicBoxes.
    let target = Target::from_coords(&water_positions(), &water_radii(), 1)
        .with_restraint(InsideBoxRestraint::new([0.0; 3], [30.0; 3], [true; 3]));
    let result = Molpack::new()
        .with_seed(42)
        .with_periodic_box([0.0; 3], [40.0; 3])
        .pack_with_report(&[target], 2);
    assert!(
        matches!(result, Err(PackError::ConflictingPeriodicBoxes { .. })),
        "expected ConflictingPeriodicBoxes, got: {result:?}"
    );
}

#[test]
fn with_periodic_box_rejects_zero_extent() {
    let (coords, radii) = single_atom();
    let target = Target::from_coords(&coords, &radii, 1);
    let result = Molpack::new()
        .with_seed(42)
        .with_periodic_box([0.0; 3], [10.0, 0.0, 10.0])
        .pack_with_report(&[target], 2);
    assert!(
        matches!(result, Err(PackError::InvalidPBCBox { .. })),
        "expected InvalidPBCBox, got: {result:?}"
    );
}

/// PBC with a shifted origin (`pbc X0 Y0 Z0  X1 Y1 Z1` script form,
/// `with_periodic_box(min, max)` API form). Verifies that:
///
/// 1. The packer accepts asymmetric, non-origin-anchored PBC bounds.
/// 2. Every packed atom lands inside the shifted reference image.
/// 3. Pairwise distances respect the minimum-image convention — i.e.
///    no two inter-molecule atom pairs are closer than `tolerance`
///    once we wrap their separation through the PBC.
///
/// Before this test, the 6-value `pbc` script form (commit 19e150f) had
/// no end-to-end coverage at the *numerical* layer — only parser-level
/// arity checks. A regression in `with_periodic_box`'s min-image wrap
/// inside the pair kernel would silently let atoms overlap across the
/// wrap boundary without tripping any test.
#[test]
fn pbc_shifted_origin_box_packs_within_bounds() {
    let min = [-10.0, -15.0, -20.0];
    let max = [20.0, 15.0, 10.0];
    let target = Target::from_coords(&water_positions(), &water_radii(), 4)
        .with_restraint(InsideBoxRestraint::new(min, max, [true; 3]));
    let tolerance: F = 2.0;
    let result = Molpack::new()
        .with_tolerance(tolerance)
        .with_seed(0xCAFE)
        .with_periodic_box(min, max)
        .pack_with_report(&[target], 5)
        .expect("shifted-origin PBC pack should succeed");

    assert_eq!(result.natoms(), 12, "expected 4 waters × 3 atoms");

    let positions = result.positions();
    let extent = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
    let half = [extent[0] * 0.5, extent[1] * 0.5, extent[2] * 0.5];

    // (1) Every atom inside the shifted reference image (allow a
    // small tolerance — the soft-confining penalty on non-periodic
    // axes leaves things a bit inside, but periodic axes wrap so
    // anything past the bound is in another image).
    let slack: F = 1e-3;
    for (i, p) in positions.iter().enumerate() {
        for k in 0..3 {
            assert!(
                p[k] >= min[k] - slack && p[k] <= max[k] + slack,
                "atom {i} axis {k} = {} out of [{}, {}]",
                p[k],
                min[k],
                max[k]
            );
        }
    }

    // (2) Pairwise minimum-image distances on inter-molecule pairs
    // (atoms 0..3 = mol 0, 3..6 = mol 1, …) must respect the
    // tolerance. We allow a small slack so the optimizer's
    // `precision` setting can leave a residual.
    let mol_of = |a: usize| a / 3;
    let mut min_dist: F = F::INFINITY;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            if mol_of(i) == mol_of(j) {
                continue;
            }
            let mut dx = positions[i][0] - positions[j][0];
            let mut dy = positions[i][1] - positions[j][1];
            let mut dz = positions[i][2] - positions[j][2];
            // Minimum-image wrap on each axis.
            if dx.abs() > half[0] {
                dx -= dx.signum() * extent[0];
            }
            if dy.abs() > half[1] {
                dy -= dy.signum() * extent[1];
            }
            if dz.abs() > half[2] {
                dz -= dz.signum() * extent[2];
            }
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            if d < min_dist {
                min_dist = d;
            }
        }
    }
    // 1e-2 matches the default Molpack `precision` floor for
    // distance-violation reporting.
    assert!(
        min_dist >= tolerance - 1e-2,
        "min-image inter-mol distance {min_dist} < tolerance {tolerance}"
    );
}

// ── error cases ────────────────────────────────────────────────────────────

#[test]
fn empty_targets_returns_error() {
    let result = Molpack::new().with_seed(42).pack_with_report(&[], 5);
    assert!(
        matches!(result, Err(PackError::NoTargets)),
        "expected NoTargets, got: {result:?}"
    );
}

// ── handler integration ────────────────────────────────────────────────────

#[test]
fn null_handler_accepted() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 2).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let result = Molpack::new()
        .with_handler(NullHandler)
        .with_seed(42)
        .pack_with_report(&[target], 5);
    assert!(result.is_ok());
}

// ── deterministic seed ─────────────────────────────────────────────────────

#[test]
fn same_seed_same_result() {
    let make_target = || {
        Target::from_coords(&water_positions(), &water_radii(), 3).with_restraint(
            InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
        )
    };
    let r1 = Molpack::new()
        .with_seed(42)
        .pack_with_report(&[make_target()], 5)
        .unwrap();
    let r2 = Molpack::new()
        .with_seed(42)
        .pack_with_report(&[make_target()], 5)
        .unwrap();
    for (a, b) in r1.positions().iter().zip(r2.positions().iter()) {
        assert!((a[0] - b[0]).abs() < 1e-6);
        assert!((a[1] - b[1]).abs() < 1e-6);
        assert!((a[2] - b[2]).abs() < 1e-6);
    }
}

// ── builder chain ──────────────────────────────────────────────────────────

#[test]
fn builder_precision_and_tolerance() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 2).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let result = Molpack::new()
        .with_precision(0.1)
        .with_tolerance(3.0)
        .with_inner_iterations(10)
        .with_seed(42)
        .pack_with_report(&[target], 5);
    assert!(result.is_ok());
}

// ── composite restraint (multiple independent restraints) ──────────────────

#[test]
fn pack_with_composite_restraints() {
    let (coords, radii) = single_atom();
    let target = Target::from_coords(&coords, &radii, 3)
        .with_restraint(InsideBoxRestraint::new(
            [-20.0, -20.0, -20.0],
            [20.0, 20.0, 20.0],
            [false; 3],
        ))
        .with_restraint(OutsideSphereRestraint::new([0.0, 0.0, 0.0], 2.0));
    let result = Molpack::new().with_seed(42).pack_with_report(&[target], 10);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 3);
}

// ── B.3: Molpack::add_restraint broadcast ──────────────────────────────────

#[test]
fn molpack_add_restraint_broadcasts_to_every_target() {
    // Spec §4 "Scope 等价律": Molpack::add_restraint(r)
    // ≡ for each target: target.with_restraint(r.clone())
    let (coords, radii) = single_atom();
    let t1 = Target::from_coords(&coords, &radii, 2).with_name("A");
    let t2 = Target::from_coords(&coords, &radii, 2).with_name("B");

    // Global restraint via Molpack — equivalent to per-target broadcast
    let result = Molpack::new()
        .with_global_restraint(InsideBoxRestraint::new([0.0; 3], [30.0; 3], [false; 3]))
        .with_seed(42)
        .pack_with_report(&[t1, t2], 5)
        .expect("global restraint should pack successfully");
    assert_eq!(result.natoms(), 4);

    // All atoms must land inside the box (all constrained)
    for pos in result.positions() {
        for k in 0..3 {
            assert!(
                pos[k] >= -0.1 && pos[k] <= 30.1,
                "global restraint should confine all targets: atom {pos:?} axis {k}"
            );
        }
    }
}

#[test]
fn molpack_add_restraint_idempotent_with_with_restraint() {
    // Both paths must produce byte-identical positions given the same seed.
    let (coords, radii) = single_atom();
    let box_r = InsideBoxRestraint::new([0.0; 3], [20.0; 3], [false; 3]);

    // Path A: global via Molpack::add_restraint
    let ta = Target::from_coords(&coords, &radii, 3);
    let ra = Molpack::new()
        .with_global_restraint(box_r)
        .with_seed(7)
        .pack_with_report(&[ta], 5)
        .unwrap();

    // Path B: per-target via with_restraint
    let tb = Target::from_coords(&coords, &radii, 3).with_restraint(box_r);
    let rb = Molpack::new()
        .with_seed(7)
        .pack_with_report(&[tb], 5)
        .unwrap();

    for (a, b) in ra.positions().iter().zip(rb.positions().iter()) {
        assert!((a[0] - b[0]).abs() < 1e-12);
        assert!((a[1] - b[1]).abs() < 1e-12);
        assert!((a[2] - b[2]).abs() < 1e-12);
    }
}

// ── avoid_overlap: initial placement around a fixed solute ───────────────────

/// Counts total GENCAN outer loops across all phases via a shared counter.
struct LoopCounter(std::sync::Arc<std::sync::atomic::AtomicUsize>);
impl molpack::Handler for LoopCounter {
    fn on_step(&mut self, _i: &molpack::StepInfo, _sys: &molpack::PackContext) {
        self.0.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }
}

/// Solvent packed around a large fixed solute converges with far less work when
/// `avoid_overlap` (default, Packmol-faithful) keeps the solvent out of the
/// solute at the initial guess. Disabling it seeds solvent atoms *inside* the
/// solute, inflating the initial overlap and forcing many extra outer loops —
/// the `pack_solvprotein` slowdown in miniature. This locks in both that the
/// `with_avoid_overlap` switch is honored and that the default reduces work.
#[test]
fn avoid_overlap_reduces_work_around_fixed_solute() {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    // Fixed solute: a dense 5×5×5 cluster (radius 1) filling the box centre,
    // whose ±1-cell exclusion covers a large fraction of the free region.
    let mut solute = Vec::new();
    for ix in -2..=2 {
        for iy in -2..=2 {
            for iz in -2..=2 {
                solute.push([ix as F * 2.0, iy as F * 2.0, iz as F * 2.0]);
            }
        }
    }
    let solute_radii = vec![1.0; solute.len()];

    let run = |avoid: bool, counter: Arc<AtomicUsize>| {
        let (c, r) = single_atom();
        let free = Target::from_coords(&c, &r, 100).with_restraint(InsideBoxRestraint::new(
            [-8.0, -8.0, -8.0],
            [8.0, 8.0, 8.0],
            [false; 3],
        ));
        let fixed = Target::from_coords(&solute, &solute_radii, 1).fixed_at([0.0, 0.0, 0.0]);
        Molpack::new()
            .with_seed(1234567)
            .with_avoid_overlap(avoid)
            .with_handler(LoopCounter(counter))
            .pack_with_report(&[free, fixed], 200)
    };

    let (on, off) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let res_on = run(true, on.clone()).expect("avoid_overlap on should converge");
    let res_off = run(false, off.clone()).expect("avoid_overlap off should converge");

    // Same valid packing either way: 100 free atoms + 125 solute atoms.
    assert_eq!(res_on.natoms(), 100 + 125);
    assert_eq!(res_off.natoms(), 100 + 125);

    // The fix's whole point: avoidance on (default) needs strictly fewer loops.
    let (loops_on, loops_off) = (on.load(Ordering::Relaxed), off.load(Ordering::Relaxed));
    assert!(
        loops_on < loops_off,
        "avoid_overlap should cut outer loops, got on={loops_on} off={loops_off}"
    );
}
