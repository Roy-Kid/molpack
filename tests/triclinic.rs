//! The `cell` builder / script surface.
//!
//! Only the declaration paths are covered here. A declared cell fixes the
//! lattice and the cell partition, but it does **not** confine anything: in
//! molpack, as in Packmol, molecules are held inside a region by a restraint,
//! and the only regions available today are axis-aligned. Until a fractional
//! `InsideCell` region exists (spec task 3), there is no way to say "inside
//! this hexagonal cell", so a pack declared with `with_cell` scatters molecules
//! across many lattice images and any minimum-image assertion over them passes
//! vacuously. Those assertions belong with task 3, against a pack that is
//! actually confined.

use molpack::{F, Molpack, Target};

const HEX_LENGTHS: [F; 3] = [26.0, 26.0, 26.0];
const HEX_ANGLES: [F; 3] = [90.0, 90.0, 120.0];

#[test]
fn a_declared_cell_and_a_periodic_box_are_mutually_exclusive() {
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 4);
    let err = Molpack::new()
        .with_cell(HEX_LENGTHS, HEX_ANGLES, [true; 3])
        .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
        .pack_with_report(&[target], 2)
        .expect_err("declaring both a cell and a box must be rejected");
    assert!(
        format!("{err}").contains("mutually exclusive"),
        "unexpected error: {err}"
    );
}

#[test]
fn a_degenerate_cell_is_rejected() {
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 4);
    // Angles that sum past the triangle inequality describe no cell.
    let err = Molpack::new()
        .with_cell([10.0, 10.0, 10.0], [170.0, 170.0, 170.0], [true; 3])
        .pack_with_report(&[target], 2)
        .expect_err("degenerate cell must be rejected");
    assert!(
        format!("{err}").contains("invalid packing cell"),
        "unexpected error: {err}"
    );
}
