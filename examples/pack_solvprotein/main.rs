//! Packmol solvprotein example: one fixed solute + water + ions in a sphere.
//!
//! Equivalent to Packmol's `solvprotein.inp`:
//! ```text
//! tolerance 2.0
//! structure protein.pdb
//!   number 1
//!   fixed 0. 0. 0. 0. 0. 0.
//!   centerofmass
//! end structure
//! structure water.pdb
//!   number 16500
//!   inside sphere 0. 0. 0. 50.
//! end structure
//! structure chloride.pdb
//!   number 20
//!   inside sphere 0. 0. 0. 50.
//! end structure
//! structure sodium.pdb
//!   number 30
//!   inside sphere 0. 0. 0. 50.
//! end structure
//! ```
//!
//! Run with:
//! ```sh
//! cargo run --release --example pack_solvprotein --features io
//! ```

use std::fs::create_dir_all;
use std::path::PathBuf;

use molpack::{
    CenteringMode, GencanPack, PackEngine, ProgressCallback, RegionRestraint, Target,
    XyzTrajectoryCallback,
};
use molrs::op::F;
use std::sync::Arc;

use molrs::core::Sphere;
use molrs::io::read_pdb;
use ndarray::array;

fn inside_sphere(center: [F; 3], radius: F) -> RegionRestraint {
    RegionRestraint(Arc::new(Sphere::new(
        array![center[0], center[1], center[2]],
        radius,
    )))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::try_init();
    let base = PathBuf::from(file!())
        .parent()
        .expect("file path has no parent")
        .to_path_buf();
    let protein = read_pdb(base.join("protein.pdb"))?;
    let water = read_pdb(base.join("water.pdb"))?;
    let sodium = read_pdb(base.join("sodium.pdb"))?;
    let chloride = read_pdb(base.join("chloride.pdb"))?;

    let sphere = inside_sphere([0.0, 0.0, 0.0], 50.0);

    let protein_target = Target::new(protein, 1)
        .with_name("protein")
        .with_centering(CenteringMode::Center)
        .fixed_at([0.0, 0.0, 0.0]);

    let water_target = Target::new(water, 16500)
        .with_restraint(sphere.clone())
        .with_name("water");

    let sodium_target = Target::new(sodium, 30)
        .with_restraint(sphere.clone())
        .with_name("sodium");

    let chloride_target = Target::new(chloride, 20)
        .with_restraint(sphere)
        .with_name("chloride");

    let mut packer = GencanPack::new();
    if std::env::var_os("MOLPACK_EXAMPLE_PROGRESS").is_some() {
        packer = packer.with_callback(Box::new(ProgressCallback::new()));
    }
    if std::env::var_os("MOLPACK_EXAMPLE_XYZ").is_some() {
        let out_dir = base.join("out");
        create_dir_all(&out_dir)?;
        packer = packer.with_callback(Box::new(XyzTrajectoryCallback::new(
            out_dir.join("solvprotein.xyz"),
            10,
        )));
    }

    let targets = vec![protein_target, water_target, sodium_target, chloride_target];
    packer.run(&targets, 800)?;

    Ok(())
}
