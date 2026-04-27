//! Lower a parsed [`Script`] to either a frame-loader-agnostic
//! [`ScriptPlan`] (no I/O) or — when the `io` feature is on — a fully
//! built [`BuildResult`] with templates already read from disk.
//!
//! Front-ends pick whichever fits:
//!
//! - **Native CLI / examples** — call [`Script::build`] (feature `io`),
//!   which reads files via molrs-io.
//! - **PyO3 / WASM / embedding hosts** — call [`Script::lower`], drive
//!   their own frame loader (e.g. molrs's Python bindings), construct
//!   [`Target`]s externally, then apply per-structure restraints with
//!   [`StructurePlan::apply`].

use std::path::{Path, PathBuf};

use crate::{
    AbovePlaneRestraint, Angle, BelowPlaneRestraint, CenteringMode, InsideBoxRestraint,
    InsideSphereRestraint, Molpack, OutsideSphereRestraint, Target,
};

use super::error::ScriptError;
use super::parser::{AtomGroup, RestraintSpec, Script, Structure};

/// Loader-agnostic lowering of a [`Script`]: paths resolved, restraints
/// kept as parsed [`RestraintSpec`]s, no filesystem access yet.
///
/// Iterate [`structures`](ScriptPlan::structures), load each
/// [`StructurePlan::filepath`] with whatever frame loader you have, build
/// a [`Target`] from the frame, and call [`StructurePlan::apply`] to
/// stamp on the script's restraints / centering / fixed placement.
pub struct ScriptPlan {
    /// Packer pre-configured with `tolerance`, `seed`, and (optional) `pbc`.
    pub packer: Molpack,
    /// One entry per `structure … end structure` block, in source order.
    pub structures: Vec<StructurePlan>,
    /// Resolved output file path.
    pub output: PathBuf,
    /// Outer-loop iteration cap (`nloop` keyword; default 400).
    pub nloop: usize,
    /// Global `filetype` override (script-level), if any. Frame loaders
    /// should fall back to extension-based detection when this is `None`.
    pub filetype: Option<String>,
}

/// Per-structure lowering: resolved file path plus the restraints / pose
/// hints that apply to its [`Target`].
pub struct StructurePlan {
    /// Absolute path to the template molecule file (resolved against the
    /// script's base directory).
    pub filepath: PathBuf,
    /// Number of copies to pack.
    pub number: usize,
    /// Molecule-wide restraints.
    pub mol_restraints: Vec<RestraintSpec>,
    /// Atom-subset restraints (`atoms … end atoms` blocks). Indices are
    /// **1-based** as written in the script; [`StructurePlan::apply`]
    /// converts to 0-based when stamping them on a [`Target`].
    pub atom_groups: Vec<AtomGroup>,
    /// Whether the `center` keyword was present.
    pub center: bool,
    /// Fixed placement: `(position [x,y,z], euler [ex,ey,ez])`.
    pub fixed: Option<([f64; 3], [f64; 3])>,
}

impl Script {
    /// Resolve paths and clone restraints into a [`ScriptPlan`] without
    /// touching the filesystem.
    ///
    /// Use this from front-ends that supply their own frame loader. The
    /// native counterpart that *does* read files is [`Script::build`]
    /// (gated behind the `io` feature).
    pub fn lower(&self, base_dir: &Path) -> Result<ScriptPlan, ScriptError> {
        if self.structures.is_empty() {
            return Err(ScriptError::NoStructures);
        }

        let mut packer = Molpack::new().with_tolerance(self.tolerance);
        if let Some(seed) = self.seed {
            packer = packer.with_seed(seed);
        }
        if let Some(pbc) = self.pbc {
            packer = packer.with_periodic_box(pbc.min, pbc.max);
        }

        let structures: Vec<StructurePlan> = self
            .structures
            .iter()
            .map(|s| StructurePlan::from_structure(s, base_dir))
            .collect();

        Ok(ScriptPlan {
            packer,
            structures,
            output: resolve(base_dir, &self.output),
            nloop: self.nloop,
            filetype: self.filetype.clone(),
        })
    }
}

impl StructurePlan {
    fn from_structure(s: &Structure, base_dir: &Path) -> Self {
        Self {
            filepath: resolve(base_dir, &s.filepath),
            number: s.number,
            mol_restraints: s.mol_restraints.clone(),
            atom_groups: s.atom_groups.clone(),
            center: s.center,
            fixed: s.fixed,
        }
    }

    /// Apply this plan's restraints / centering / fixed pose onto a
    /// [`Target`] the caller built from the template frame.
    pub fn apply(&self, mut target: Target) -> Target {
        for r in &self.mol_restraints {
            target = apply_mol_restraint(target, r);
        }

        for group in &self.atom_groups {
            target = apply_atom_group(target, group);
        }

        if self.center {
            target = target.with_centering(CenteringMode::Center);
        }

        if let Some((pos, euler)) = self.fixed {
            target = target.fixed_at(pos).with_orientation([
                Angle::from_radians(euler[0]),
                Angle::from_radians(euler[1]),
                Angle::from_radians(euler[2]),
            ]);
        }

        target
    }
}

fn resolve(base: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
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

// ─────────────────────────────────────────────────────────────────────
// `io`-gated convenience: load every structure file via molrs-io.
// ─────────────────────────────────────────────────────────────────────

/// Everything a script expanded to: a configured packer, the target
/// list, the resolved output path, and the outer-loop iteration cap.
///
/// The packer is not yet equipped with a handler; callers decide
/// whether to attach a [`ProgressHandler`](crate::ProgressHandler),
/// a custom handler, or none.
#[cfg(feature = "io")]
pub struct BuildResult {
    pub packer: Molpack,
    pub targets: Vec<Target>,
    pub output: PathBuf,
    pub nloop: usize,
}

#[cfg(feature = "io")]
impl Script {
    /// Lower the script *and* read each template via molrs-io.
    ///
    /// Available when the `io` feature is on; equivalent to calling
    /// [`Script::lower`] then loading each structure file with
    /// [`super::io::read_frame`] and applying its [`StructurePlan`].
    pub fn build(&self, base_dir: &Path) -> Result<BuildResult, ScriptError> {
        let plan = self.lower(base_dir)?;
        let filetype = plan.filetype.as_deref();
        let targets: Vec<Target> = plan
            .structures
            .iter()
            .map(|sp| -> Result<Target, ScriptError> {
                let frame = super::io::read_frame(&sp.filepath, filetype)?;
                Ok(sp.apply(Target::new(frame, sp.number)))
            })
            .collect::<Result<_, _>>()?;
        Ok(BuildResult {
            packer: plan.packer,
            targets,
            output: plan.output,
            nloop: plan.nloop,
        })
    }
}
