//! Packmol mixture example: water + urea in a box.
//!
//! Equivalent to Packmol's `mixture.inp`:
//! ```text
//! tolerance 2.0
//! output mixture.pdb
//! structure water.pdb
//!   number 1000
//!   inside box 0. 0. 0. 40. 40. 40.
//! end structure
//! structure urea.pdb
//!   number 400
//!   inside box 0. 0. 0. 40. 40. 40.
//! end structure
//! ```
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_mixture --features io
//! ```

use std::fs::create_dir_all;
use std::path::PathBuf;

use molpack::{
    GencanPack, PackEngine, ProgressCallback, RegionRestraint, Target, XyzTrajectoryCallback,
};
use std::sync::Arc;

use molrs::core::Cuboid;
use molrs::io::read_pdb;
use ndarray::array;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();
    let base = PathBuf::from(file!())
        .parent()
        .expect("file path has no parent")
        .to_path_buf();
    let water = read_pdb(base.join("water.pdb"))?;
    let urea = read_pdb(base.join("urea.pdb"))?;

    let box_restraint = RegionRestraint(Arc::new(Cuboid::new(
        array![0.0, 0.0, 0.0],
        array![40.0, 40.0, 40.0],
    )));

    let water_target = Target::new(water, 1000)
        .with_restraint(box_restraint.clone())
        .with_name("water");

    let urea_target = Target::new(urea, 400)
        .with_restraint(box_restraint)
        .with_name("urea");

    let mut packer = GencanPack::new();
    if std::env::var_os("MOLPACK_EXAMPLE_PROGRESS").is_some() {
        packer = packer.with_callback(Box::new(ProgressCallback::new()));
    }
    if std::env::var_os("MOLPACK_EXAMPLE_XYZ").is_some() {
        let out_dir = base.join("out");
        create_dir_all(&out_dir)?;
        packer = packer.with_callback(Box::new(XyzTrajectoryCallback::new(
            out_dir.join("mixture.xyz"),
            10,
        )));
    }

    let targets = vec![water_target, urea_target];
    packer.run(&targets, 400)?;

    Ok(())
}
