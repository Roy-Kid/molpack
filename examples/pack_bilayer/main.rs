//! Packmol bilayer example: a double layer with solvent above and below.
//!
//! Based on Packmol's `bilayer.inp` from https://m3g.github.io/packmol/examples.shtml.
//! This version maps Packmol atom-level orientation constraints directly:
//! - atoms 31 32 below plane z=2
//! - atoms 1 2 over plane z=12
//! - atoms 1 2 below plane z=16
//! - atoms 31 32 over plane z=26
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_bilayer --features io
//! ```

use std::fs::create_dir_all;
use std::path::PathBuf;

use molpack::{GenCanPack, PackEngine, ProgressHandler, RegionRestraint, Target, XYZHandler};
use molrs::op::types::F;
use std::sync::Arc;

use molrs::io::data::pdb::read_pdb_frame;
use molrs::spatial::region::{Cuboid, HalfSpace, NotRegion};
use ndarray::array;

// ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

fn inside_box(min: [F; 3], max: [F; 3]) -> RegionRestraint {
    RegionRestraint(Arc::new(Cuboid::new(
        array![min[0], min[1], min[2]],
        array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
    )))
}

/// `n · x >= d`: the complement of the half-space behind the plane.
fn above_plane(normal: [F; 3], distance: F) -> RegionRestraint {
    let n = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let point = [
        distance * normal[0] / n,
        distance * normal[1] / n,
        distance * normal[2] / n,
    ];
    RegionRestraint(Arc::new(NotRegion::new(Arc::new(
        HalfSpace::new(normal, point).expect("plane"),
    ))))
}

/// `n · x <= d`: the half-space behind the plane.
fn below_plane(normal: [F; 3], distance: F) -> RegionRestraint {
    let n = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let point = [
        distance * normal[0] / n,
        distance * normal[1] / n,
        distance * normal[2] / n,
    ];
    RegionRestraint(Arc::new(HalfSpace::new(normal, point).expect("plane")))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();
    let base = PathBuf::from(file!())
        .parent()
        .expect("file path has no parent")
        .to_path_buf();
    let water = read_pdb_frame(base.join("water.pdb"))?;
    let lipid = read_pdb_frame(base.join("palmitoil.pdb"))?;

    let water_low = Target::new(water.clone(), 50)
        .with_restraint(inside_box([0.0, 0.0, -10.0], [40.0, 40.0, 0.0]))
        .with_name("water_low");

    let water_high = Target::new(water, 50)
        .with_restraint(inside_box([0.0, 0.0, 28.0], [40.0, 40.0, 38.0]))
        .with_name("water_high");

    let lipid_low = Target::new(lipid.clone(), 10)
        .with_restraint(inside_box([0.0, 0.0, 0.0], [40.0, 40.0, 14.0]))
        .with_atom_restraint(&[30, 31], below_plane([0.0, 0.0, 1.0], 2.0))
        .with_atom_restraint(&[0, 1], above_plane([0.0, 0.0, 1.0], 12.0))
        .with_name("lipid_low");

    let lipid_high = Target::new(lipid, 10)
        .with_restraint(inside_box([0.0, 0.0, 14.0], [40.0, 40.0, 28.0]))
        .with_atom_restraint(&[0, 1], below_plane([0.0, 0.0, 1.0], 16.0))
        .with_atom_restraint(&[30, 31], above_plane([0.0, 0.0, 1.0], 26.0))
        .with_name("lipid_high");

    let targets = vec![water_low, water_high, lipid_low, lipid_high];
    let mut packer = GenCanPack::new();
    if std::env::var_os("MOLPACK_EXAMPLE_PROGRESS").is_some() {
        packer = packer.with_handler(Box::new(ProgressHandler::new()));
    }
    if std::env::var_os("MOLPACK_EXAMPLE_XYZ").is_some() {
        let out_dir = base.join("out");
        create_dir_all(&out_dir)?;
        packer = packer.with_handler(Box::new(XYZHandler::new(out_dir.join("bilayer.xyz"), 10)));
    }

    packer.run(&targets, 800)?;

    Ok(())
}
