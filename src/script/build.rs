//! Lower a parsed [`Script`] to a ready-to-run [`Molpack`] + [`Target`] set.

use std::path::{Path, PathBuf};

use crate::{
    AbovePlaneRestraint, Angle, BelowPlaneRestraint, CenteringMode, InsideBoxRestraint,
    InsideSphereRestraint, Molpack, OutsideSphereRestraint, Target,
};

use super::error::ScriptError;
use super::io::read_frame;
use super::parser::{AtomGroup, RestraintSpec, Script, Structure};

/// Everything a script expanded to: a configured packer, the target
/// list, the resolved output path, and the outer-loop iteration cap.
///
/// The packer is not yet equipped with a handler; callers decide
/// whether to attach a [`ProgressHandler`](crate::ProgressHandler),
/// a custom handler, or none.
pub struct BuildResult {
    pub packer: Molpack,
    pub targets: Vec<Target>,
    pub output: PathBuf,
    pub nloop: usize,
}

impl Script {
    /// Convert the script into a configured packer and target list.
    ///
    /// Relative paths inside the script (structure files, output) are
    /// resolved against `base_dir`.
    pub fn build(&self, base_dir: &Path) -> Result<BuildResult, ScriptError> {
        if self.structures.is_empty() {
            return Err(ScriptError::NoStructures);
        }

        let mut packer = Molpack::new().with_tolerance(self.tolerance);
        if let Some(seed) = self.seed {
            packer = packer.with_seed(seed);
        }

        let targets: Vec<Target> = self
            .structures
            .iter()
            .map(|s| build_target(s, self.filetype.as_deref(), base_dir))
            .collect::<Result<_, _>>()?;

        Ok(BuildResult {
            packer,
            targets,
            output: resolve(base_dir, &self.output),
            nloop: self.nloop,
        })
    }
}

fn resolve(base: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

fn build_target(
    s: &Structure,
    global_filetype: Option<&str>,
    base_dir: &Path,
) -> Result<Target, ScriptError> {
    let path = resolve(base_dir, &s.filepath);
    let frame = read_frame(&path, global_filetype)?;
    let mut target = Target::new(frame, s.number);

    for r in &s.mol_restraints {
        target = apply_mol_restraint(target, r);
    }

    for group in &s.atom_groups {
        target = apply_atom_group(target, group);
    }

    if s.center {
        target = target.with_centering(CenteringMode::Center);
    }

    if let Some((pos, euler)) = s.fixed {
        target = target.fixed_at(pos).with_orientation([
            Angle::from_radians(euler[0]),
            Angle::from_radians(euler[1]),
            Angle::from_radians(euler[2]),
        ]);
    }

    Ok(target)
}

fn apply_mol_restraint(target: Target, r: &RestraintSpec) -> Target {
    match *r {
        RestraintSpec::InsideBox { min, max } => {
            target.with_restraint(InsideBoxRestraint::new(min, max, [false; 3]))
        }
        RestraintSpec::InsideSphere { center, radius } => {
            target.with_restraint(InsideSphereRestraint::new(center, radius))
        }
        RestraintSpec::OutsideSphere { center, radius } => {
            target.with_restraint(OutsideSphereRestraint::new(center, radius))
        }
        RestraintSpec::AbovePlane { normal, distance } => {
            target.with_restraint(AbovePlaneRestraint::new(normal, distance))
        }
        RestraintSpec::BelowPlane { normal, distance } => {
            target.with_restraint(BelowPlaneRestraint::new(normal, distance))
        }
    }
}

fn apply_atom_group(mut target: Target, group: &AtomGroup) -> Target {
    // Script indices are 1-based; Target::with_atom_restraint expects 0-based.
    let zero_indexed: Vec<usize> = group
        .atom_indices
        .iter()
        .map(|&i| i.saturating_sub(1))
        .collect();
    let indices = zero_indexed.as_slice();
    for r in &group.restraints {
        target = match *r {
            RestraintSpec::InsideBox { min, max } => {
                target.with_atom_restraint(indices, InsideBoxRestraint::new(min, max, [false; 3]))
            }
            RestraintSpec::InsideSphere { center, radius } => {
                target.with_atom_restraint(indices, InsideSphereRestraint::new(center, radius))
            }
            RestraintSpec::OutsideSphere { center, radius } => {
                target.with_atom_restraint(indices, OutsideSphereRestraint::new(center, radius))
            }
            RestraintSpec::AbovePlane { normal, distance } => {
                target.with_atom_restraint(indices, AbovePlaneRestraint::new(normal, distance))
            }
            RestraintSpec::BelowPlane { normal, distance } => {
                target.with_atom_restraint(indices, BelowPlaneRestraint::new(normal, distance))
            }
        };
    }
    target
}
