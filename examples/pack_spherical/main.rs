//! Packmol spherical example: double-layered shell with water inside and outside.
//!
//! Based on Packmol's `spherical.inp` from https://m3g.github.io/packmol/examples.shtml.
//!
//! Packmol input structure (full counts: 308 / 90 / 300 / 17536):
//! ```text
//! structure water.pdb              # inner water
//!   number 308
//!   inside sphere 0. 0. 0. 13.
//! end structure
//!
//! structure palmitoil.pdb           # inner layer
//!   number 90
//!   atoms 37
//!     inside sphere 0. 0. 0. 14.
//!   end atoms
//!   atoms 5
//!     outside sphere 0. 0. 0. 26.
//!   end atoms
//! end structure
//!
//! structure palmitoil.pdb           # outer layer
//!   number 300
//!   atoms 5
//!     inside sphere 0. 0. 0. 29.
//!   end atoms
//!   atoms 37
//!     outside sphere 0. 0. 0. 41.
//!   end atoms
//! end structure
//!
//! structure water.pdb              # outer water
//!   number 17536
//!   inside box -47.5 -47.5 -47.5 47.5 47.5 47.5
//!   outside sphere 0. 0. 0. 43.
//! end structure
//! ```
//!
//! Layer molecules use `with_restraint`, same semantics as all other restraints.
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_spherical --features io
//! ```

use std::path::PathBuf;

use molpack::{F, GenCanPack, PackEngine, ProgressHandler, RegionRestraint, Target};
use std::sync::Arc;

use molrs::io::data::pdb::read_pdb_frame;
use molrs::spatial::region::{Cuboid, NotRegion, Sphere};
use ndarray::array;

// ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

fn inside_box(min: [F; 3], max: [F; 3]) -> RegionRestraint {
    RegionRestraint(Arc::new(Cuboid::new(
        array![min[0], min[1], min[2]],
        array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
    )))
}

fn inside_sphere(center: [F; 3], radius: F) -> RegionRestraint {
    RegionRestraint(Arc::new(Sphere::new(
        array![center[0], center[1], center[2]],
        radius,
    )))
}

fn outside_sphere(center: [F; 3], radius: F) -> RegionRestraint {
    RegionRestraint(Arc::new(NotRegion::new(Arc::new(Sphere::new(
        array![center[0], center[1], center[2]],
        radius,
    )))))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();
    let base = PathBuf::from(file!())
        .parent()
        .expect("file path has no parent")
        .to_path_buf();
    let water = read_pdb_frame(base.join("water.pdb"))?;
    let lipid = read_pdb_frame(base.join("palmitoil.pdb"))?;

    let origin = [0.0, 0.0, 0.0];

    // 1. Inner water sphere: 308 molecules inside sphere r=13
    let water_inner = Target::new(water.clone(), 308)
        .with_restraint(inside_sphere(origin, 13.0))
        .with_name("water_inner");

    // 2. Inner layer: 90 molecules.
    //    Packmol input constrains only specific atoms:
    //    atom 37 inside sphere r=14, atom 5 outside sphere r=26.
    let lipid_inner = Target::new(lipid.clone(), 90)
        .with_atom_restraint(&[36], inside_sphere(origin, 14.0))
        .with_atom_restraint(&[4], outside_sphere(origin, 26.0))
        .with_name("lipid_inner");

    // 3. Outer layer: 300 molecules.
    //    Packmol input constrains only specific atoms:
    //    atom 5 inside sphere r=29, atom 37 outside sphere r=41.
    let lipid_outer = Target::new(lipid, 300)
        .with_atom_restraint(&[4], inside_sphere(origin, 29.0))
        .with_atom_restraint(&[36], outside_sphere(origin, 41.0))
        .with_name("lipid_outer");

    // 4. Outer water shell: 17536 molecules, box ±47.5, outside sphere r=43
    let water_outer = Target::new(water, 17536)
        .with_restraint(inside_box([-47.5, -47.5, -47.5], [47.5, 47.5, 47.5]))
        .with_restraint(outside_sphere(origin, 43.0))
        .with_name("water_outer");

    // Target order matches Packmol: water_inner → lipid_inner → lipid_outer → water_outer
    let targets = vec![water_inner, lipid_inner, lipid_outer, water_outer];
    let mut packer = GenCanPack::new();
    if std::env::var_os("MOLPACK_EXAMPLE_PROGRESS").is_some() {
        packer = packer.with_handler(Box::new(ProgressHandler::new()));
    }

    // Match spherical-comment.inp defaults:
    // - nloop defaults to 200 * ntype (ntype = 4 => 800)
    // - seed defaults to 1234567
    packer.run(&targets, 800)?;

    Ok(())
}
