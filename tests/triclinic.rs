//! Packing into a non-orthorhombic cell.
//!
//! Packmol supports orthorhombic periodic boundaries only; a hexagonal or
//! triclinic cell has to be approximated by an enclosing box and cropped, or
//! expanded into a rectangular supercell. These tests cover the geometry that
//! approximation stands in for.
//!
//! The oracle is deliberately hand-rolled: minimum distances are recomputed by
//! stepping over the 27 lattice images directly from `a`, `b`, `c`, sharing no
//! code with the packer's cell list. A check routed through the same partition
//! would be wrong in exactly the cases the partition is wrong.

use molpack::{
    AbovePlaneRestraint, BelowPlaneRestraint, F, GenCanPack, InsideCellRegion, PackEngine,
    RegionRestraint, Target,
};

/// Lattice vectors (columns of H) for a cell given by lengths and angles, in
/// the same upper-triangular convention `SimBox` uses.
fn lattice(lengths: [F; 3], angles_deg: [F; 3]) -> [[F; 3]; 3] {
    let [a, b, c] = lengths;
    let [alpha, beta, gamma] = angles_deg.map(F::to_radians);
    let xy = b * gamma.cos();
    let xz = c * beta.cos();
    let ly = (b * b - xy * xy).sqrt();
    let yz = (b * c * alpha.cos() - xy * xz) / ly;
    let lz = (c * c - xz * xz - yz * yz).sqrt();
    [[a, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]]
}

/// Fractional coordinates of a Cartesian point, by back-substitution through
/// the upper-triangular lattice.
fn fractional(p: [F; 3], h: [[F; 3]; 3]) -> [F; 3] {
    let fz = p[2] / h[2][2];
    let fy = (p[1] - h[1][2] * fz) / h[1][1];
    let fx = (p[0] - h[0][1] * fy - h[0][2] * fz) / h[0][0];
    [fx, fy, fz]
}

/// Smallest distance between atoms of different molecules, minimised over the
/// lattice images of the periodic axes.
fn min_intermolecular_distance(coords: &[[F; 3]], h: [[F; 3]; 3], pbc: [bool; 3]) -> F {
    let range = |periodic: bool| if periodic { -1..=1 } else { 0..=0 };
    let mut best = F::INFINITY;
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            let d0 = [
                coords[j][0] - coords[i][0],
                coords[j][1] - coords[i][1],
                coords[j][2] - coords[i][2],
            ];
            for na in range(pbc[0]) {
                for nb in range(pbc[1]) {
                    for nc in range(pbc[2]) {
                        let shift = [na as F, nb as F, nc as F];
                        let mut d = d0;
                        for (axis, item) in d.iter_mut().enumerate() {
                            *item += h[axis][0] * shift[0]
                                + h[axis][1] * shift[1]
                                + h[axis][2] * shift[2];
                        }
                        let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                        if dist < best {
                            best = dist;
                        }
                    }
                }
            }
        }
    }
    best
}

const HEX_LENGTHS: [F; 3] = [26.0, 26.0, 26.0];
const HEX_ANGLES: [F; 3] = [90.0, 90.0, 120.0];

fn pack_in_cell(
    lengths: [F; 3],
    angles: [F; 3],
    pbc: [bool; 3],
    n: usize,
    tolerance: F,
    seed: u64,
) -> Vec<[F; 3]> {
    let cell = InsideCellRegion::from_lengths_angles(lengths, angles, pbc).expect("cell");
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[tolerance / 2.0], n)
        .with_restraint(RegionRestraint(cell));
    GenCanPack::new()
        .with_seed(seed)
        .with_tolerance(tolerance)
        .run(&[target], 30)
        .expect("pack")
        .positions()
        .to_vec()
}

#[test]
fn hexagonal_pack_stays_within_the_cell_to_within_the_tolerance() {
    // `Inside*` restraints are quadratic penalties, not hard walls, so
    // equilibrium leaves sub-tolerance excursions past a face — measured worst
    // case here is ~1.1 Å against a 2 Å tolerance. Under periodicity such an
    // excursion is not even a defect: fractional 1.04 and 0.04 are the same
    // configuration and the minimum image cannot tell them apart (which is why
    // the distance checks below hold regardless). What this test guards is that
    // molecules stay *near* the declared cell at all: with nothing confining
    // them they drift across hundreds of lattice images.
    let tolerance: F = 2.0;
    let coords = pack_in_cell(HEX_LENGTHS, HEX_ANGLES, [true; 3], 200, tolerance, 42);
    assert_eq!(coords.len(), 200);
    let h = lattice(HEX_LENGTHS, HEX_ANGLES);
    let spacing = [h[0][0], h[1][1], h[2][2]];
    for p in &coords {
        let f = fractional(*p, h);
        for k in 0..3 {
            let outside = (-f[k]).max(f[k] - 1.0).max(0.0) * spacing[k];
            assert!(
                outside <= tolerance,
                "molecule sits {outside:.3} Å past the cell face on axis {k}, \
                 beyond the {tolerance} Å tolerance"
            );
        }
    }
}

#[test]
fn hexagonal_pack_satisfies_tolerance_under_the_true_minimum_image() {
    let tolerance: F = 2.0;
    let coords = pack_in_cell(HEX_LENGTHS, HEX_ANGLES, [true; 3], 200, tolerance, 42);
    let h = lattice(HEX_LENGTHS, HEX_ANGLES);
    let dmin = min_intermolecular_distance(&coords, h, [true; 3]);
    assert!(
        dmin >= tolerance - 1e-6,
        "closest pair {dmin:.4} A is below the {tolerance} A tolerance once the tilted \
         lattice images are taken into account"
    );
    // The check has to bite: a pack loose enough that no pair comes near the
    // tolerance would satisfy the assertion above without exercising the
    // minimum image at all.
    assert!(
        dmin < 2.0 * tolerance,
        "pack is too dilute ({dmin:.4} A) for the minimum-image check to mean anything"
    );
}

#[test]
fn strongly_tilted_pack_satisfies_tolerance() {
    let lengths = [24.0, 26.0, 28.0];
    let angles = [70.0, 80.0, 65.0];
    let tolerance: F = 2.0;
    let coords = pack_in_cell(lengths, angles, [true; 3], 200, tolerance, 7);
    let h = lattice(lengths, angles);
    let dmin = min_intermolecular_distance(&coords, h, [true; 3]);
    assert!(
        dmin >= tolerance - 1e-6,
        "closest pair {dmin:.4} A is below the {tolerance} A tolerance"
    );
    assert!(dmin < 2.0 * tolerance, "pack too dilute: {dmin:.4} A");
}

#[test]
fn slab_cell_is_periodic_in_plane_and_confined_along_the_normal() {
    // xy periodic, z confined: the shape of every interface system. The oracle
    // wraps only the two periodic axes, so a pack that wrapped along z shows up
    // here as a violation.
    let lengths = [24.0, 24.0, 30.0];
    let angles = [90.0, 90.0, 120.0];
    let tolerance: F = 2.0;
    let coords = pack_in_cell(lengths, angles, [true, true, false], 150, tolerance, 11);
    let h = lattice(lengths, angles);
    let dmin = min_intermolecular_distance(&coords, h, [true, true, false]);
    assert!(
        dmin >= tolerance - 1e-6,
        "closest pair {dmin:.4} A is below the {tolerance} A tolerance"
    );
    assert!(dmin < 2.0 * tolerance, "pack too dilute: {dmin:.4} A");
}

#[test]
fn the_cell_is_declared_once_by_the_region_alone() {
    // No `with_cell` on the builder: the packer picks the lattice up from the
    // region that confines the molecules.
    let coords = pack_in_cell(HEX_LENGTHS, HEX_ANGLES, [true; 3], 60, 2.0, 5);
    assert_eq!(coords.len(), 60);
}

#[test]
fn a_declared_cell_and_a_periodic_box_are_mutually_exclusive() {
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 4);
    let err = GenCanPack::new()
        .with_cell(HEX_LENGTHS, HEX_ANGLES, [true; 3])
        .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
        .run(&[target], 2)
        .expect_err("declaring both a cell and a box must be rejected");
    assert!(
        format!("{err}").contains("mutually exclusive"),
        "unexpected error: {err}"
    );
}

#[test]
fn a_degenerate_cell_is_rejected() {
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 4);
    let err = GenCanPack::new()
        .with_cell([10.0, 10.0, 10.0], [170.0, 170.0, 170.0], [true; 3])
        .run(&[target], 2)
        .expect_err("degenerate cell must be rejected");
    assert!(
        format!("{err}").contains("invalid packing cell"),
        "unexpected error: {err}"
    );
}

// ── half-spaces under periodicity ──────────────────────────────────────────

#[test]
fn a_plane_across_a_periodic_axis_is_rejected() {
    // z is periodic, so a plane with a z-component has no well-defined side:
    // translating by c moves a point across it.
    let cell =
        InsideCellRegion::from_lengths_angles(HEX_LENGTHS, HEX_ANGLES, [true; 3]).expect("cell");
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 8)
        .with_restraint(RegionRestraint(cell))
        .with_restraint(AbovePlaneRestraint::new([0.0, 0.0, 1.0], 5.0));
    let err = GenCanPack::new()
        .with_seed(1)
        .run(&[target], 2)
        .expect_err("plane across a periodic axis must be rejected");
    let msg = format!("{err}");
    assert!(
        msg.contains("periodic lattice direction 2"),
        "the error must name the offending axis, got: {msg}"
    );
}

#[test]
fn a_plane_along_a_confined_axis_is_accepted() {
    // Same plane, but with z non-periodic: now it is exactly the slab
    // constraint an interface system needs.
    let cell = InsideCellRegion::from_lengths_angles(HEX_LENGTHS, HEX_ANGLES, [true, true, false])
        .expect("cell");
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 8)
        .with_restraint(RegionRestraint(cell))
        .with_restraint(AbovePlaneRestraint::new([0.0, 0.0, 1.0], 5.0));
    let result = GenCanPack::new().with_seed(1).run(&[target], 10);
    assert!(
        result.is_ok(),
        "expected the slab plane to be accepted: {result:?}"
    );
}

#[test]
fn a_plane_in_the_periodic_plane_is_rejected_naming_the_first_axis() {
    // A normal lying in the periodic xy plane crosses lattice vector a.
    let cell = InsideCellRegion::from_lengths_angles(HEX_LENGTHS, HEX_ANGLES, [true, true, false])
        .expect("cell");
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 8)
        .with_restraint(RegionRestraint(cell))
        .with_restraint(BelowPlaneRestraint::new([1.0, 0.0, 0.0], 5.0));
    let err = GenCanPack::new()
        .with_seed(1)
        .run(&[target], 2)
        .expect_err("plane across a periodic axis must be rejected");
    assert!(
        format!("{err}").contains("periodic lattice direction 0"),
        "unexpected error: {err}"
    );
}
