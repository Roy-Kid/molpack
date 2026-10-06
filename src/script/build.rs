//! Lower a parsed [`Script`] to either a frame-loader-agnostic
//! [`ScriptPlan`] (no I/O) or — when the `io` feature is on — a fully
//! built `BuildResult` with templates already read from disk.
//!
//! Front-ends pick whichever fits:
//!
//! - **Native CLI / examples** — call `Script::build` (feature `io`),
//!   which reads files via molrs-io.
//! - **PyO3 / WASM / embedding hosts** — call [`Script::lower`], drive
//!   their own frame loader (e.g. molrs's Python bindings), construct
//!   [`Target`]s externally, then apply per-structure restraints with
//!   [`StructurePlan::apply`].

use std::path::{Path, PathBuf};

use crate::restraint::geometric::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideBoxRestraint, InsideCubeRestraint,
    InsideCylinderRestraint, InsideEllipsoidRestraint, InsideSphereRestraint, OutsideBoxRestraint,
    OutsideCubeRestraint, OutsideCylinderRestraint, OutsideEllipsoidRestraint,
    OutsideSphereRestraint,
};
use crate::{Angle, AtomRestraint, CenteringMode, GenCanPack, PackEngine, Target};

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
    /// Engine pre-configured with `tolerance`, `seed`, and (optional) `pbc`.
    pub entry: GenCanPack,
    /// One entry per `structure … end structure` block, in source order.
    pub structures: Vec<StructurePlan>,
    /// Resolved output file path.
    pub output: PathBuf,
    /// Outer-loop iteration cap (`nloop` keyword; default `200 * ntype`).
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
    /// Structure-level `radius`, applied to every atom before any
    /// atom-specific override.
    pub radius: Option<f64>,
    /// Structure-level `fscale`.
    pub fscale: Option<f64>,
    /// Structure-level `short_radius`.
    pub short_radius: Option<f64>,
    /// Structure-level `short_radius_scale`.
    pub short_radius_scale: Option<f64>,
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
    /// native counterpart that *does* read files is `Script::build`,
    /// which is compiled only with the `io` feature on — hence a plain
    /// code span here instead of a cross-reference, since a
    /// default-feature documentation build has no such item to link to.
    pub fn lower(&self, base_dir: &Path) -> Result<ScriptPlan, ScriptError> {
        if self.structures.is_empty() {
            return Err(ScriptError::NoStructures);
        }

        let mut entry = GenCanPack::new()
            .with_tolerance(self.tolerance)
            .with_avoid_overlap(self.avoid_overlap);
        if let Some(seed) = self.seed {
            entry = entry.with_seed(seed);
        }
        if let Some(pbc) = self.pbc {
            entry = entry.with_periodic_box(pbc.min, pbc.max, [true; 3]);
        }
        if let Some(cell) = self.cell {
            entry = entry.with_cell(cell.lengths, cell.angles_deg, cell.pbc);
        }

        let structures: Vec<StructurePlan> = self
            .structures
            .iter()
            .map(|s| StructurePlan::from_structure(s, base_dir))
            .collect();

        Ok(ScriptPlan {
            entry,
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
            radius: s.radius,
            fscale: s.fscale,
            short_radius: s.short_radius,
            short_radius_scale: s.short_radius_scale,
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

        // Radii before restraints of the same scope, and structure level
        // before atom level — Packmol runs the two passes in that order
        // (`app/packmol.f90` lines 294 and 390) so the narrower selection wins.
        if let Some(r) = self.radius {
            target = target.with_radius(r);
        }
        if let Some(v) = self.fscale {
            target = target.with_fscale(v);
        }
        if let Some(v) = self.short_radius {
            target = target.with_short_radius(v);
        }
        if let Some(v) = self.short_radius_scale {
            target = target.with_short_radius_scale(v);
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

/// Lower one parsed [`RestraintSpec`] to its concrete [`AtomRestraint`].
///
/// Single source of truth for the spec → restraint mapping, shared by the
/// whole-molecule and per-atom-group application paths. Adding a new restraint
/// kind is one arm here plus one parser arm and one [`RestraintSpec`] variant.
fn restraint_from_spec(r: &RestraintSpec) -> Box<dyn AtomRestraint> {
    match *r {
        RestraintSpec::InsideBox { min, max } => Box::new(InsideBoxRestraint::new(min, max)),
        RestraintSpec::OutsideBox { min, max } => Box::new(OutsideBoxRestraint::new(min, max)),
        RestraintSpec::InsideCube { origin, side } => {
            Box::new(InsideCubeRestraint::new(origin, side))
        }
        RestraintSpec::OutsideCube { origin, side } => {
            Box::new(OutsideCubeRestraint::new(origin, side))
        }
        RestraintSpec::InsideSphere { center, radius } => {
            Box::new(InsideSphereRestraint::new(center, radius))
        }
        RestraintSpec::OutsideSphere { center, radius } => {
            Box::new(OutsideSphereRestraint::new(center, radius))
        }
        RestraintSpec::InsideEllipsoid {
            center,
            axes,
            exponent,
        } => Box::new(InsideEllipsoidRestraint::new(center, axes, exponent)),
        RestraintSpec::OutsideEllipsoid {
            center,
            axes,
            exponent,
        } => Box::new(OutsideEllipsoidRestraint::new(center, axes, exponent)),
        RestraintSpec::InsideCylinder {
            center,
            axis,
            radius,
            length,
        } => Box::new(InsideCylinderRestraint::new(center, axis, radius, length)),
        RestraintSpec::OutsideCylinder {
            center,
            axis,
            radius,
            length,
        } => Box::new(OutsideCylinderRestraint::new(center, axis, radius, length)),
        RestraintSpec::AbovePlane { normal, distance } => {
            Box::new(AbovePlaneRestraint::new(normal, distance))
        }
        RestraintSpec::BelowPlane { normal, distance } => {
            Box::new(BelowPlaneRestraint::new(normal, distance))
        }
    }
}

fn apply_mol_restraint(target: Target, r: &RestraintSpec) -> Target {
    target.with_restraint(restraint_from_spec(r))
}

fn apply_atom_group(mut target: Target, group: &AtomGroup) -> Target {
    // Script indices are 1-based; Target::with_atom_restraint expects 0-based.
    let zero_indexed: Vec<usize> = group
        .atom_indices
        .iter()
        .map(|&i| i.saturating_sub(1))
        .collect();
    let indices = zero_indexed.as_slice();
    if let Some(r) = group.radius {
        target = target.with_atom_radius(indices, r);
    }
    if let Some(v) = group.fscale {
        target = target.with_atom_fscale(indices, v);
    }
    if let Some(v) = group.short_radius {
        target = target.with_atom_short_radius(indices, v);
    }
    if let Some(v) = group.short_radius_scale {
        target = target.with_atom_short_radius_scale(indices, v);
    }
    for r in &group.restraints {
        target = target.with_atom_restraint(indices, restraint_from_spec(r));
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
    pub entry: GenCanPack,
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
    /// `molrs::io::read_frame` (the script's `filetype`, else the file
    /// name, picks the format) and applying its [`StructurePlan`].
    pub fn build(&self, base_dir: &Path) -> Result<BuildResult, ScriptError> {
        let plan = self.lower(base_dir)?;
        let filetype = plan.filetype.as_deref();
        let targets: Vec<Target> = plan
            .structures
            .iter()
            .map(|sp| -> Result<Target, ScriptError> {
                let frame =
                    molrs::io::read_frame(&sp.filepath, filetype).map_err(|e| ScriptError::Io {
                        path: sp.filepath.clone(),
                        message: format!("reading template: {e}"),
                    })?;
                Ok(sp.apply(Target::new(frame, sp.number)))
            })
            .collect::<Result<_, _>>()?;
        Ok(BuildResult {
            entry: plan.entry,
            targets,
            output: plan.output,
            nloop: plan.nloop,
        })
    }
}

#[cfg(test)]
mod atom_property_tests {
    //! The four per-atom packing properties at both `.inp` levels —
    //! structure keyword and `atoms ... end atoms` block — plus the
    //! lowering that layers them onto a [`Target`](crate::Target).

    use crate::Target;
    use crate::script::parse;

    fn plan(src: &str) -> crate::script::ScriptPlan {
        parse(src)
            .expect("parse")
            .lower(std::path::Path::new("."))
            .expect("lower")
    }

    /// Structure-level `radius`, as Packmol reads it outside an `atoms` block.
    #[test]
    fn structure_radius_is_parsed() {
        let p = plan(
            "tolerance 2.0\noutput o.xyz\n\
                 structure a.pdb\n  number 3\n  radius 3.5\n\
                 inside box 0. 0. 0. 10. 10. 10.\nend structure\n",
        );
        assert_eq!(p.structures[0].radius, Some(3.5));
    }

    /// Atom-specific `radius`, inside an `atoms ... end atoms` block.
    #[test]
    fn atom_group_radius_is_parsed() {
        let p = plan(
            "tolerance 2.0\noutput o.xyz\n\
                 structure a.pdb\n  number 1\n\
                 inside box 0. 0. 0. 10. 10. 10.\n\
                 atoms 1 3\n  radius 6.0\nend atoms\n\
                 end structure\n",
        );
        let g = &p.structures[0].atom_groups[0];
        assert_eq!(g.atom_indices, vec![1, 3]);
        assert_eq!(g.radius, Some(6.0));
    }

    /// Lowering must reproduce the Rust API's layering, with the script's
    /// 1-based indices mapped to 0-based.
    #[test]
    fn lowering_applies_both_radius_layers() {
        let p = plan(
            "tolerance 4.0\noutput o.xyz\n\
                 structure a.pdb\n  number 2\n  radius 3.0\n\
                 inside box 0. 0. 0. 10. 10. 10.\n\
                 atoms 3\n  radius 6.0\nend atoms\n\
                 end structure\n",
        );
        let t = p.structures[0].apply(crate::Target::from_coords(
            &[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.5; 3],
            2,
        ));
        assert_eq!(t.resolved_radii(2.0), vec![3.0, 3.0, 6.0]);
    }

    #[test]
    fn a_non_positive_radius_is_a_parse_error() {
        let src = "tolerance 2.0\noutput o.xyz\n\
                       structure a.pdb\n  number 1\n  radius -1.0\n\
                       inside box 0. 0. 0. 10. 10. 10.\nend structure\n";
        assert!(parse(src).is_err(), "negative radius must be rejected");
    }
    fn plan_body(body: &str) -> crate::script::ScriptPlan {
        let src = format!(
            "tolerance 4.0\noutput o.xyz\nstructure a.pdb\n  number 1\n\
                 inside box 0. 0. 0. 10. 10. 10.\n{body}end structure\n"
        );
        parse(&src)
            .expect("parse")
            .lower(std::path::Path::new("."))
            .expect("lower")
    }

    #[test]
    fn structure_level_keywords_are_parsed() {
        let p = plan_body("  fscale 0.5\n  short_radius 0.75\n  short_radius_scale 4.0\n");
        let s = &p.structures[0];
        assert_eq!(s.fscale, Some(0.5));
        assert_eq!(s.short_radius, Some(0.75));
        assert_eq!(s.short_radius_scale, Some(4.0));
    }

    #[test]
    fn atom_level_keywords_are_parsed() {
        let p = plan_body("atoms 2\n  fscale 0.5\n  short_radius 0.75\nend atoms\n");
        let g = &p.structures[0].atom_groups[0];
        assert_eq!(g.fscale, Some(0.5));
        assert_eq!(g.short_radius, Some(0.75));
    }

    #[test]
    fn lowering_applies_every_layer() {
        let p = plan_body(
            "  fscale 0.5\n  short_radius 0.75\n\
                 atoms 3\n  fscale 2.0\n  short_radius_scale 9.0\nend atoms\n",
        );
        let t = p.structures[0].apply(Target::from_coords(
            &[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.5; 3],
            1,
        ));
        assert_eq!(t.resolved_fscale(), vec![0.5, 0.5, 2.0]);
        assert_eq!(t.resolved_short_radii(1.0), vec![0.75; 3]);
        assert_eq!(t.resolved_short_radius_scale(3.0), vec![3.0, 3.0, 9.0]);
        assert_eq!(t.uses_short_radius(), vec![true; 3]);
    }

    #[test]
    fn a_non_positive_fscale_is_a_parse_error() {
        let src = "tolerance 4.0\noutput o.xyz\nstructure a.pdb\n  number 1\n  fscale 0\n\
                       inside box 0. 0. 0. 10. 10. 10.\nend structure\n";
        assert!(parse(src).is_err());
    }
}
