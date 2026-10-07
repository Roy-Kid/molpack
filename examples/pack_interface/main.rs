//! Packmol interface example: water/chloroform interface with fixed molecule.
//!
//! Equivalent to Packmol's `interface.inp`:
//! ```text
//! tolerance 2.0
//! output interface.xyz
//! structure water.xyz
//!   number 100
//!   inside box -20. 0. 0. 0. 39. 39.
//! end structure
//! structure chlor.xyz
//!   number 30
//!   inside box 0. 0. 0. 21. 39. 39.
//! end structure
//! structure t3.xyz
//!   center
//!   fixed 0. 20. 20. 1.57 1.57 1.57
//! end structure
//! ```
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_interface --features io
//! ```

use std::fs::create_dir_all;
use std::path::PathBuf;

use molpack::{
    Angle, CenteringMode, GenCanPack, PackEngine, ProgressHandler, RegionRestraint, Target,
    XYZHandler,
};
use molrs::op::F;
use std::sync::Arc;

use molrs::core::Cuboid;
use molrs::io::read_pdb;
use ndarray::array;

// ── molrs regions lifted to "stay inside" (the one geometric restraint) ─────

fn inside_box(min: [F; 3], max: [F; 3]) -> RegionRestraint {
    RegionRestraint(Arc::new(Cuboid::new(
        array![min[0], min[1], min[2]],
        array![max[0] - min[0], max[1] - min[1], max[2] - min[2]],
    )))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();
    let base = PathBuf::from(file!())
        .parent()
        .expect("file path has no parent")
        .to_path_buf();
    let water = read_pdb(base.join("water.pdb"))?;
    let chloroform = read_pdb(base.join("chloroform.pdb"))?;
    let t3 = read_pdb(base.join("t3.pdb"))?;

    let water_target = Target::new(water, 100)
        .with_restraint(inside_box([-20.0, 0.0, 0.0], [0.0, 39.0, 39.0]))
        .with_name("water");

    let chloro_target = Target::new(chloroform, 30)
        .with_restraint(inside_box([0.0, 0.0, 0.0], [21.0, 39.0, 39.0]))
        .with_name("chloroform");

    let t3_target = Target::new(t3, 1)
        .with_name("t3")
        .with_centering(CenteringMode::Center)
        .fixed_at([0.0, 20.0, 20.0])
        .with_orientation([
            Angle::from_radians(1.57),
            Angle::from_radians(1.57),
            Angle::from_radians(1.57),
        ]);

    let mut packer = GenCanPack::new();
    if std::env::var_os("MOLPACK_EXAMPLE_PROGRESS").is_some() {
        packer = packer.with_handler(Box::new(ProgressHandler::new()));
    }
    if std::env::var_os("MOLPACK_EXAMPLE_XYZ").is_some() {
        let out_dir = base.join("out");
        create_dir_all(&out_dir)?;
        packer = packer.with_handler(Box::new(XYZHandler::new(out_dir.join("interface.xyz"), 10)));
    }

    let targets = vec![water_target, chloro_target, t3_target];
    packer.run(&targets, 400)?;

    Ok(())
}
