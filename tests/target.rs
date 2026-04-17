//! Tests for Target builder: construction, natoms/count, fixed_at,
//! centering modes, restraint attachment, and hook validation.

use molpack::{F, InsideBoxRestraint, InsideSphereRestraint, Molpack, Target};

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

// ── construction ───────────────────────────────────────────────────────────

#[test]
fn from_coords_basic() {
    let t = Target::from_coords(&water_positions(), &water_radii(), 5);
    assert_eq!(t.natoms(), 3);
    assert_eq!(t.count, 5);
    assert!(t.name.is_none());
    assert!(t.fixed_at.is_none());
}

#[test]
fn with_name() {
    let t = Target::from_coords(&water_positions(), &water_radii(), 1).with_name("water");
    assert_eq!(t.name, Some("water".to_string()));
}

#[test]
fn ref_coords_are_centered() {
    let coords = vec![[10.0, 20.0, 30.0], [12.0, 20.0, 30.0]];
    let t = Target::from_coords(&coords, &[1.0, 1.0], 1);
    // Center should be at (11, 20, 30), so ref_coords are [-1, 0, 0] and [1, 0, 0]
    let cx: F = t.ref_coords.iter().map(|p| p[0]).sum::<F>() / t.natoms() as F;
    let cy: F = t.ref_coords.iter().map(|p| p[1]).sum::<F>() / t.natoms() as F;
    let cz: F = t.ref_coords.iter().map(|p| p[2]).sum::<F>() / t.natoms() as F;
    assert!(cx.abs() < 1e-6, "ref_coords x center should be 0");
    assert!(cy.abs() < 1e-6, "ref_coords y center should be 0");
    assert!(cz.abs() < 1e-6, "ref_coords z center should be 0");
}

#[test]
fn input_coords_preserved() {
    let coords = vec![[10.0, 20.0, 30.0]];
    let t = Target::from_coords(&coords, &[1.0], 1);
    assert!((t.input_coords[0][0] - 10.0).abs() < 1e-6);
    assert!((t.input_coords[0][1] - 20.0).abs() < 1e-6);
}

#[test]
fn new_uses_geometric_center_even_when_elements_are_known() {
    use molrs::block::Block;
    use ndarray::Array1;

    let mut atoms = Block::new();
    atoms
        .insert("x", Array1::from_vec(vec![0.0, 1.0]).into_dyn())
        .expect("insert x");
    atoms
        .insert("y", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
        .expect("insert y");
    atoms
        .insert("z", Array1::from_vec(vec![0.0, 0.0]).into_dyn())
        .expect("insert z");
    atoms
        .insert(
            "element",
            Array1::from_vec(vec!["O".to_string(), "H".to_string()]).into_dyn(),
        )
        .expect("insert element");

    let mut frame = molrs::Frame::new();
    frame.insert("atoms", atoms);

    let t = Target::new(frame, 1);
    let arithmetic_center = (t.ref_coords[0][0] + t.ref_coords[1][0]) / 2.0;
    assert!(
        arithmetic_center.abs() < 1e-6,
        "geometry center should be zero"
    );
}

// ── restraints ─────────────────────────────────────────────────────────────

#[test]
fn with_restraint() {
    let t = Target::from_coords(&water_positions(), &water_radii(), 5).with_restraint(
        InsideBoxRestraint::new([0.0, 0.0, 0.0], [20.0, 20.0, 20.0], [false; 3]),
    );
    assert_eq!(t.molecule_restraints.len(), 1);
}

#[test]
fn with_restraint_chained() {
    let t = Target::from_coords(&water_positions(), &water_radii(), 5)
        .with_restraint(InsideBoxRestraint::new(
            [0.0, 0.0, 0.0],
            [20.0, 20.0, 20.0],
            [false; 3],
        ))
        .with_restraint(InsideSphereRestraint::new([10.0, 10.0, 10.0], 50.0));
    assert_eq!(t.molecule_restraints.len(), 2);
}

#[test]
fn with_atom_restraint() {
    // Indices are now 0-based (matching Rust convention) — no internal
    // conversion happens. Caller subtracts 1 when porting from Packmol
    // `.inp` files.
    let t = Target::from_coords(&water_positions(), &water_radii(), 5)
        .with_atom_restraint(&[0, 1], InsideSphereRestraint::new([0.0, 0.0, 0.0], 5.0));
    assert_eq!(t.atom_restraints.len(), 1);
    assert_eq!(t.atom_restraints[0].0, vec![0, 1]);
}

// ── fixed placement ────────────────────────────────────────────────────────

#[test]
fn fixed_at_sets_count_to_1() {
    // Explicit count=1 to satisfy the assertion in fixed_at().
    let t = Target::from_coords(&water_positions(), &water_radii(), 1).fixed_at([0.0, 0.0, 0.0]);
    assert_eq!(t.count, 1);
    assert!(t.fixed_at.is_some());
    let fp = t.fixed_at.unwrap();
    assert!((fp.position[0]).abs() < 1e-6);
    assert!((fp.orientation[0].radians()).abs() < 1e-6);
}

#[test]
fn fixed_at_with_orientation() {
    use molpack::Angle;
    let t = Target::from_coords(&water_positions(), &water_radii(), 1)
        .fixed_at([1.0, 2.0, 3.0])
        .with_orientation([
            Angle::from_radians(0.1),
            Angle::from_radians(0.2),
            Angle::from_radians(0.3),
        ]);
    assert_eq!(t.count, 1);
    let fp = t.fixed_at.unwrap();
    assert!((fp.position[0] - 1.0).abs() < 1e-6);
    assert!((fp.orientation[2].radians() - 0.3).abs() < 1e-6);
}

#[test]
fn fixed_target_auto_centering_disabled() {
    // When fixed_at is used with Auto centering (default), the fixed molecule
    // should NOT be centered — its input coords are used directly.
    let free = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1).with_restraint(
        InsideBoxRestraint::new([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], [false; 3]),
    );
    let fixed = Target::from_coords(&[[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]], &[1.0, 1.0], 1)
        .fixed_at([0.0, 0.0, 0.0]);

    let result = Molpack::new()
        .with_seed(1)
        .pack(&[free, fixed], 5)
        .expect("pack should succeed");

    // Fixed atoms follow free atoms in output.
    assert!((result.positions()[1][0] - 10.0).abs() < 1e-6);
    assert!((result.positions()[2][0] - 12.0).abs() < 1e-6);
}

#[test]
fn fixed_target_centered() {
    let free = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1).with_restraint(
        InsideBoxRestraint::new([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0], [false; 3]),
    );
    let fixed = Target::from_coords(&[[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]], &[1.0, 1.0], 1)
        .with_centering(molpack::CenteringMode::Center)
        .fixed_at([0.0, 0.0, 0.0]);

    let result = Molpack::new()
        .with_seed(1)
        .pack(&[free, fixed], 5)
        .expect("pack should succeed");

    // COM of [10,12] = 11. After centering, ref_coords = [-1, +1].
    // Placed at origin → positions = [-1, +1].
    assert!((result.positions()[1][0] + 1.0).abs() < 1e-6);
    assert!((result.positions()[2][0] - 1.0).abs() < 1e-6);
}

// ── centering modes ────────────────────────────────────────────────────────

#[test]
fn with_centering_center() {
    use molpack::CenteringMode;
    let t = Target::from_coords(&water_positions(), &water_radii(), 1)
        .with_centering(CenteringMode::Center);
    assert_eq!(t.centering, CenteringMode::Center);
}

#[test]
fn with_centering_off() {
    use molpack::CenteringMode;
    let t = Target::from_coords(&water_positions(), &water_radii(), 1)
        .with_centering(CenteringMode::Off);
    assert_eq!(t.centering, CenteringMode::Off);
}

// ── rotation constraints ───────────────────────────────────────────────────

#[test]
fn with_rotation_bound() {
    use molpack::{Angle, Axis};
    let t = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1)
        .with_rotation_bound(Axis::X, Angle::from_degrees(0.0), Angle::from_degrees(10.0))
        .with_rotation_bound(Axis::Y, Angle::from_degrees(90.0), Angle::from_degrees(5.0))
        .with_rotation_bound(
            Axis::Z,
            Angle::from_degrees(180.0),
            Angle::from_degrees(15.0),
        );
    // Euler variable order: [beta(Y) = 0, gama(Z) = 1, teta(X) = 2]
    assert!(t.rotation_bound[0].is_some()); // beta (Y)
    assert!(t.rotation_bound[1].is_some()); // gama (Z)
    assert!(t.rotation_bound[2].is_some()); // teta (X)
}

// ── perturb budget ─────────────────────────────────────────────────────────

#[test]
fn with_perturb_budget() {
    let t = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 1).with_perturb_budget(5);
    assert_eq!(t.perturb_budget, Some(5));
}

// ── default element ────────────────────────────────────────────────────────

#[test]
fn default_elements_are_x() {
    let t = Target::from_coords(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &[1.0, 1.0], 1);
    assert_eq!(t.elements, vec!["X", "X"]);
}

// ── panics ─────────────────────────────────────────────────────────────────

#[test]
#[should_panic(expected = "positions and radii must have the same length")]
fn mismatched_coords_and_radii_panics() {
    Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0, 2.0], 1);
}
