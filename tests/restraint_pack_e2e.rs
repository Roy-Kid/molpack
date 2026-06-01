//! End-to-end pack tests for the 9 restraint kinds that the gated
//! `examples_batch.rs` Packmol-fixture suite never exercises:
//!
//! | Restraint                      | Packmol kind |
//! |--------------------------------|-------------:|
//! | `InsideCubeRestraint`          |            2 |
//! | `InsideEllipsoidRestraint`     |            5 |
//! | `OutsideCubeRestraint`         |            6 |
//! | `OutsideBoxRestraint`          |            7 |
//! | `OutsideEllipsoidRestraint`    |            9 |
//! | `InsideCylinderRestraint`      |           12 |
//! | `OutsideCylinderRestraint`     |           13 |
//! | `AboveGaussianRestraint`       |           14 |
//! | `BelowGaussianRestraint`       |           15 |
//!
//! Each test packs a handful of single-atom molecules through the full
//! `Molpack::pack` driver and then runs the corresponding analytic
//! restraint over the result to assert that every packed atom lies in
//! the satisfied region. `validate_from_targets` also catches pairwise
//! tolerance violations.
//!
//! Speed budget: ≤ 200 ms per test (debug build); tiny atom counts and a
//! capped `max_loops` keep this in the default tier.

use molpack::{
    AboveGaussianRestraint, BelowGaussianRestraint, F, InsideBoxRestraint, InsideCubeRestraint,
    InsideCylinderRestraint, InsideEllipsoidRestraint, Molpack, OutsideBoxRestraint,
    OutsideCubeRestraint, OutsideCylinderRestraint, OutsideEllipsoidRestraint, Restraint, Target,
    validate_from_targets,
};

const SEED: u64 = 0xC0FFEE;
const TOL: F = 2.0;
const PRECISION: F = 1e-2;

// ── helpers ────────────────────────────────────────────────────────────────

fn single_atom() -> (Vec<[F; 3]>, Vec<F>) {
    (vec![[0.0, 0.0, 0.0]], vec![1.0])
}

/// Run `Molpack::pack` and validate atom-count + distance tolerance.
fn pack_and_validate(targets: Vec<Target>, max_loops: usize) {
    let result = Molpack::new()
        .with_tolerance(TOL)
        .with_precision(PRECISION)
        .with_seed(SEED)
        .pack(&targets, max_loops)
        .expect("pack should succeed");
    let positions = result.positions().to_vec();
    let report = validate_from_targets(&targets, &positions, TOL, PRECISION);
    assert!(
        report.is_valid(),
        "validate_from_targets failed: {report:?}"
    );
}

/// Assert that every packed atom satisfies a single-restraint analytic
/// check (`f(pos) ≤ slack²`). Run AFTER `pack_and_validate` so atom
/// count is already verified.
fn assert_all_atoms_satisfy<R: Restraint>(positions: &[[F; 3]], r: &R, label: &str) {
    let slack: F = 1.0; // slacker than `precision` — restraints decay smoothly.
    for (i, p) in positions.iter().enumerate() {
        let v = r.f(p, 1.0, 0.01);
        assert!(v <= slack, "{label} atom {i} = {p:?} unsatisfied: f={v}");
    }
}

/// Pack one target and return its final positions for restraint
/// satisfaction checking, while also running `validate_from_targets`.
fn pack_one(target: Target, max_loops: usize) -> Vec<[F; 3]> {
    let result = Molpack::new()
        .with_tolerance(TOL)
        .with_precision(PRECISION)
        .with_seed(SEED)
        .pack(std::slice::from_ref(&target), max_loops)
        .expect("pack should succeed");
    let positions = result.positions().to_vec();
    let report = validate_from_targets(&[target], &positions, TOL, PRECISION);
    assert!(
        report.is_valid(),
        "validate_from_targets failed: {report:?}"
    );
    positions
}

// ── kind 2: InsideCubeRestraint ────────────────────────────────────────────

#[test]
fn pack_inside_cube_e2e() {
    let (coords, radii) = single_atom();
    let cube = InsideCubeRestraint::new([0.0, 0.0, 0.0], 20.0);
    let target = Target::from_coords(&coords, &radii, 5).with_restraint(cube);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &cube, "inside_cube");
}

// ── kind 5: InsideEllipsoidRestraint ───────────────────────────────────────

#[test]
fn pack_inside_ellipsoid_e2e() {
    let (coords, radii) = single_atom();
    let ell = InsideEllipsoidRestraint::new([0.0, 0.0, 0.0], [12.0, 8.0, 6.0], 1.0);
    let target = Target::from_coords(&coords, &radii, 4).with_restraint(ell);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &ell, "inside_ellipsoid");
}

// ── kind 6: OutsideCubeRestraint (paired with confining outer box) ─────────

#[test]
fn pack_outside_cube_e2e() {
    let (coords, radii) = single_atom();
    // Atoms must stay inside outer 30 Å box AND outside inner 8 Å cube.
    let outer = InsideBoxRestraint::new([0.0; 3], [30.0; 3], [false; 3]);
    let inner = OutsideCubeRestraint::new([10.0, 10.0, 10.0], 10.0);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(inner);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &inner, "outside_cube");
}

// ── kind 7: OutsideBoxRestraint ────────────────────────────────────────────

#[test]
fn pack_outside_box_e2e() {
    let (coords, radii) = single_atom();
    let outer = InsideBoxRestraint::new([0.0; 3], [30.0; 3], [false; 3]);
    let inner = OutsideBoxRestraint::new([10.0, 10.0, 10.0], [20.0, 20.0, 20.0]);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(inner);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &inner, "outside_inner_box");
}

// ── kind 9: OutsideEllipsoidRestraint ──────────────────────────────────────

/// Note: kind 9 has a known `f`-vs-`fg` `scale2` asymmetry (`f` ignores
/// `scale2`, gradient applies it). End-to-end pack still converges
/// because the optimizer only sees the gradient; this test locks in
/// pack-level satisfaction without touching that internal detail.
#[test]
fn pack_outside_ellipsoid_e2e() {
    let (coords, radii) = single_atom();
    let outer = InsideBoxRestraint::new([-15.0; 3], [15.0; 3], [false; 3]);
    let inner = OutsideEllipsoidRestraint::new([0.0, 0.0, 0.0], [4.0, 3.0, 2.5], 1.0);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(inner);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &inner, "outside_ellipsoid");
}

// ── kind 12: InsideCylinderRestraint ───────────────────────────────────────

#[test]
fn pack_inside_cylinder_e2e() {
    let (coords, radii) = single_atom();
    // Cylinder along +z, base at origin, length 20, radius 6.
    let cyl = InsideCylinderRestraint::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 6.0, 20.0);
    let target = Target::from_coords(&coords, &radii, 4).with_restraint(cyl);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &cyl, "inside_cylinder");
}

// ── kind 13: OutsideCylinderRestraint ──────────────────────────────────────

#[test]
fn pack_outside_cylinder_e2e() {
    let (coords, radii) = single_atom();
    let outer = InsideBoxRestraint::new([-15.0; 3], [15.0; 3], [false; 3]);
    // Inner cylinder along z, radius 3, length 20.
    let cyl = OutsideCylinderRestraint::new([0.0, 0.0, -10.0], [0.0, 0.0, 1.0], 3.0, 20.0);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(cyl);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &cyl, "outside_cylinder");
}

// ── kind 14: AboveGaussianRestraint ────────────────────────────────────────

#[test]
fn pack_above_gaussian_e2e() {
    let (coords, radii) = single_atom();
    let outer = InsideBoxRestraint::new([-15.0, -15.0, 0.0], [15.0, 15.0, 30.0], [false; 3]);
    // Gaussian centred at origin, σ=4, base z₀=0, height 5.
    let gauss = AboveGaussianRestraint::new(0.0, 0.0, 4.0, 4.0, 0.0, 5.0);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(gauss);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &gauss, "above_gaussian");
}

// ── kind 15: BelowGaussianRestraint ────────────────────────────────────────

#[test]
fn pack_below_gaussian_e2e() {
    let (coords, radii) = single_atom();
    let outer = InsideBoxRestraint::new([-15.0, -15.0, -30.0], [15.0, 15.0, 0.0], [false; 3]);
    // Gaussian centred at origin, σ=4, base z₀=0, height -5 (mirror).
    let gauss = BelowGaussianRestraint::new(0.0, 0.0, 4.0, 4.0, 0.0, -5.0);
    let target = Target::from_coords(&coords, &radii, 4)
        .with_restraint(outer)
        .with_restraint(gauss);

    let positions = pack_one(target, 5);
    assert_all_atoms_satisfy(&positions, &outer, "inside_outer_box");
    assert_all_atoms_satisfy(&positions, &gauss, "below_gaussian");
}

// ── multi-target sanity (cube + ellipsoid in one box) ──────────────────────

/// Drive `pack_and_validate` (multi-target validator path) on two
/// distinct restraint kinds that the gated suite never sees together.
/// This is the only test in this file that uses the multi-target
/// `validate_from_targets` helper; everything else goes through
/// `pack_one`.
#[test]
fn pack_cube_and_ellipsoid_multi_target_e2e() {
    let (coords, radii) = single_atom();
    let t1 = Target::from_coords(&coords, &radii, 3)
        .with_name("in_cube")
        .with_restraint(InsideCubeRestraint::new([0.0, 0.0, 0.0], 20.0));
    let t2 = Target::from_coords(&coords, &radii, 3)
        .with_name("in_ellipsoid")
        .with_restraint(InsideEllipsoidRestraint::new(
            [10.0, 10.0, 10.0],
            [8.0, 6.0, 5.0],
            1.0,
        ));
    pack_and_validate(vec![t1, t2], 5);
}
