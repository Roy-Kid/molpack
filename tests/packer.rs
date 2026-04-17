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
    let result = Molpack::new().with_seed(42).pack(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 3);
}

#[test]
fn pack_three_waters_in_box() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 3).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [40.0, 40.0, 40.0], [false; 3]),
    );
    let result = Molpack::new().with_seed(42).pack(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 9);
}

#[test]
fn pack_in_sphere() {
    let target = Target::from_coords(&water_positions(), &water_radii(), 3)
        .with_restraint(InsideSphereRestraint::new([0.0, 0.0, 0.0], 20.0));
    let result = Molpack::new().with_seed(42).pack(&[target], 5);
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
    let result = Molpack::new().with_seed(42).pack(&[t1, t2], 5);
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
        .pack(&[free, fixed], 5)
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
    let result = Molpack::new().with_seed(42).pack(&[target], 5);
    assert!(result.is_ok());
    assert_eq!(result.unwrap().natoms(), 6);
}

#[test]
fn invalid_pbc_box_rejected() {
    // Zero-extent y-axis on a periodic restraint must be rejected.
    let target = Target::from_coords(&water_positions(), &water_radii(), 1).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [10.0, 0.0, 10.0], [true; 3]),
    );
    let result = Molpack::new().with_seed(7).pack(&[target], 5);
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
    let result = Molpack::new().with_seed(42).pack(&[t1, t2], 5);
    assert!(
        matches!(result, Err(PackError::ConflictingPeriodicBoxes { .. })),
        "expected ConflictingPeriodicBoxes, got: {result:?}"
    );
}

// ── error cases ────────────────────────────────────────────────────────────

#[test]
fn empty_targets_returns_error() {
    let result = Molpack::new().with_seed(42).pack(&[], 5);
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
        .pack(&[target], 5);
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
        .pack(&[make_target()], 5)
        .unwrap();
    let r2 = Molpack::new()
        .with_seed(42)
        .pack(&[make_target()], 5)
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
        .pack(&[target], 5);
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
    let result = Molpack::new().with_seed(42).pack(&[target], 10);
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
        .pack(&[t1, t2], 5)
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
        .pack(&[ta], 5)
        .unwrap();

    // Path B: per-target via with_restraint
    let tb = Target::from_coords(&coords, &radii, 3).with_restraint(box_r);
    let rb = Molpack::new().with_seed(7).pack(&[tb], 5).unwrap();

    for (a, b) in ra.positions().iter().zip(rb.positions().iter()) {
        assert!((a[0] - b[0]).abs() < 1e-12);
        assert!((a[1] - b[1]).abs() < 1e-12);
        assert!((a[2] - b[2]).abs() < 1e-12);
    }
}
