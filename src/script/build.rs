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

use crate::restraint::profile::{
    Coordinate, DensityFloor, Distribution, InputKind, ProfilePenalty, ProfileRestraint,
    ProfileTarget, ShellJacobian, TabulatedProfile,
};
use crate::{
    AbovePlaneRestraint, Angle, BelowPlaneRestraint, CenteringMode, InsideBoxRestraint,
    InsideCubeRestraint, InsideCylinderRestraint, InsideEllipsoidRestraint, InsideSphereRestraint,
    Molpack, OutsideBoxRestraint, OutsideCubeRestraint, OutsideCylinderRestraint,
    OutsideEllipsoidRestraint, OutsideSphereRestraint, Restraint, Target,
};

use super::error::ScriptError;
use super::parser::{AtomGroup, RestraintSpec, Script, Structure};
use super::parser_profile::{ProfileDistribution, ProfileGeometry, ProfileInputKind};

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

        let mut packer = Molpack::new()
            .with_tolerance(self.tolerance)
            .with_avoid_overlap(self.avoid_overlap);
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

/// Lower one parsed [`RestraintSpec`] to its concrete [`Restraint`].
///
/// Single source of truth for the spec → restraint mapping, shared by the
/// whole-molecule and per-atom-group application paths. Adding a new restraint
/// kind is one arm here plus one parser arm and one [`RestraintSpec`] variant.
fn restraint_from_spec(r: &RestraintSpec) -> Box<dyn Restraint> {
    match *r {
        RestraintSpec::InsideBox { min, max } => {
            Box::new(InsideBoxRestraint::new(min, max, [false; 3]))
        }
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
        RestraintSpec::Profile {
            ref dist,
            geometry,
            input_kind,
        } => Box::new(lower_profile(dist, geometry, input_kind)),
    }
}

// ─────────────────────────────────────────────────────────────────────
// `profile` keyword lowering
// ─────────────────────────────────────────────────────────────────────

/// Fixed energy-scale factor kT (shared with the overlap penalty's units) used
/// for every script-lowered profile restraint; sets the bias steepness.
const PROFILE_KT: crate::F = 1.0;
/// Inner clamp radius (Å) for radial/cylindrical coordinates, keeping ∇ξ finite
/// at the singularity.
const PROFILE_R_GUARD: crate::F = 1e-3;
/// Density floor ρ_min and reference ρ₀ for the energy cap, turning a
/// zero-density bin into a finite plateau (no Inf/NaN gradient at pack time).
const PROFILE_RHO_MIN: crate::F = 1e-6;
const PROFILE_RHO0: crate::F = 1.0;

/// Lower a parsed `profile` spec into a concrete [`ProfileRestraint`].
///
/// Total by construction: the parser has already validated every numeric input
/// (positive σ / width / λ, non-zero normal / axis, a well-formed tabulated
/// grid), so the coordinate and tabulated-profile constructors here cannot fail.
fn lower_profile(
    dist: &ProfileDistribution,
    geometry: ProfileGeometry,
    input_kind: ProfileInputKind,
) -> ProfileRestraint {
    let coordinate = lower_coordinate(geometry);
    let jacobian = lower_jacobian(geometry);
    let kind = lower_input_kind(input_kind);
    let floor = DensityFloor::new(PROFILE_RHO_MIN, PROFILE_RHO0, PROFILE_KT);
    let target = lower_target(dist, jacobian, kind, floor);
    ProfileRestraint::non_periodic(coordinate, target, PROFILE_KT)
}

/// Map the parsed geometry to a [`Coordinate`]. The parser rejects a zero
/// normal/axis, so the fallible constructors succeed.
fn lower_coordinate(geometry: ProfileGeometry) -> Coordinate {
    match geometry {
        ProfileGeometry::Plane { normal, point } => {
            Coordinate::planar(normal, point).expect("parser rejects a zero profile-plane normal")
        }
        ProfileGeometry::Radial { center } => Coordinate::radial(center, PROFILE_R_GUARD),
        ProfileGeometry::Cylinder { origin, axis, .. } => {
            Coordinate::cylindrical(origin, axis, PROFILE_R_GUARD)
                .expect("parser rejects a zero profile-cylinder axis")
        }
    }
}

/// Map the parsed geometry to its shell-volume [`ShellJacobian`].
fn lower_jacobian(geometry: ProfileGeometry) -> ShellJacobian {
    match geometry {
        ProfileGeometry::Plane { .. } => ShellJacobian::Planar,
        ProfileGeometry::Radial { .. } => ShellJacobian::Radial,
        ProfileGeometry::Cylinder { length, .. } => ShellJacobian::Cylindrical { length },
    }
}

/// Map the parser's input-kind flag to the restraint layer's [`InputKind`].
fn lower_input_kind(input_kind: ProfileInputKind) -> InputKind {
    match input_kind {
        ProfileInputKind::Density => InputKind::VolumetricDensity,
        ProfileInputKind::Histogram => InputKind::CountHistogram,
    }
}

/// Build the [`ProfileTarget`] — an analytic [`ProfilePenalty`] for the
/// closed-form shapes or a [`TabulatedProfile`] for an inline node table.
fn lower_target(
    dist: &ProfileDistribution,
    jacobian: ShellJacobian,
    input_kind: InputKind,
    floor: DensityFloor,
) -> ProfileTarget {
    match dist {
        ProfileDistribution::Gaussian { mu, sigma } => analytic(
            Distribution::Gaussian {
                mu: *mu,
                sigma: *sigma,
            },
            jacobian,
            input_kind,
            floor,
        ),
        ProfileDistribution::Erf { xi0, w, rising } => analytic(
            Distribution::ErfStep {
                xi0: *xi0,
                w: *w,
                rising: *rising,
            },
            jacobian,
            input_kind,
            floor,
        ),
        ProfileDistribution::Tanh { xi0, w, rising } => analytic(
            Distribution::TanhStep {
                xi0: *xi0,
                w: *w,
                rising: *rising,
            },
            jacobian,
            input_kind,
            floor,
        ),
        ProfileDistribution::Exponential { lambda } => analytic(
            Distribution::Exponential { lambda: *lambda },
            jacobian,
            input_kind,
            floor,
        ),
        ProfileDistribution::Tabulated { nodes } => {
            let xs: Vec<crate::F> = nodes.iter().map(|&(xi, _)| xi).collect();
            let vs: Vec<crate::F> = nodes.iter().map(|&(_, rho)| rho).collect();
            let profile = TabulatedProfile::new(xs, vs, input_kind, jacobian, floor)
                .expect("parser validates the tabulated grid before lowering");
            ProfileTarget::Tabulated(profile)
        }
    }
}

/// Compose an analytic [`Distribution`] into a boxed [`ProfileTarget::Analytic`].
fn analytic(
    dist: Distribution,
    jacobian: ShellJacobian,
    input_kind: InputKind,
    floor: DensityFloor,
) -> ProfileTarget {
    ProfileTarget::Analytic(ProfilePenalty {
        dist,
        jacobian,
        input_kind,
        floor,
    })
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::F;

    /// Parse a single `profile` line inside a structure block and lower its
    /// restraint through `restraint_from_spec`.
    fn lower_one_profile(line: &str) -> Box<dyn Restraint> {
        let src =
            format!("output out.pdb\n\nstructure mol.pdb\n  number 1\n  {line}\nend structure\n");
        let script = crate::script::parse(&src).expect("parse failed");
        let spec = script.structures[0].mol_restraints[0].clone();
        restraint_from_spec(&spec)
    }

    // ── ac-002: a lowered profile spec reproduces the analytic f / fg ────────

    #[test]
    fn lower_gaussian_plane_matches_harmonic_form() {
        // Planar coordinate n̂ = ẑ through the origin → ξ = z, so a Gaussian
        // profile is the harmonic well U = (kt/2σ²)(z−μ)² with kt = PROFILE_KT.
        let r = lower_one_profile("profile gaussian plane 0. 0. 1. 0. 0. 0. mu 10. sigma 2.");
        let (mu, sigma): (F, F) = (10.0, 2.0);
        let scale: F = 1.0;
        let k = PROFILE_KT / (sigma * sigma); // spring constant
        for &z in &[6.0, 8.0, 10.0, 13.0] {
            let x = [0.0, 0.0, z];
            let want_u = scale * 0.5 * k * (z - mu).powi(2);
            let got_u = r.f(&x, scale, 0.0);
            assert!(
                (got_u - want_u).abs() <= 1e-9,
                "U @ z={z}: got {got_u}, want {want_u}"
            );
            // Gradient: dU/dz · ẑ = scale·k·(z−μ) on the z-axis, zero on x/y.
            let mut g = [0.0; 3];
            let got_u2 = r.fg(&x, scale, 0.0, &mut g);
            assert!((got_u2 - want_u).abs() <= 1e-9, "fg value mismatch @ z={z}");
            let want_gz = scale * k * (z - mu);
            assert!((g[0]).abs() <= 1e-12 && (g[1]).abs() <= 1e-12);
            assert!(
                (g[2] - want_gz).abs() <= 1e-9,
                "dU/dz @ z={z}: got {}, want {want_gz}",
                g[2]
            );
        }
    }

    #[test]
    fn lower_exponential_radial_density_has_constant_inward_force() {
        // Radial coordinate centred at origin → ξ = r; an exponential *density*
        // (density flag → no shell correction) gives a constant |dU/dξ| = kt/λ
        // pointing outward in +ξ, i.e. along +x̂ on the +x axis.
        let r = lower_one_profile("profile exponential radial 0. 0. 0. lambda 2.5 density");
        let lambda: F = 2.5;
        let want_mag = PROFILE_KT / lambda; // 0.4
        let x = [5.0, 0.0, 0.0];
        let mut g = [0.0; 3];
        let _ = r.fg(&x, 1.0, 0.0, &mut g);
        assert!(
            (g[0] - want_mag).abs() <= 1e-9,
            "radial exp force gx: got {}, want {want_mag}",
            g[0]
        );
        assert!(g[1].abs() <= 1e-12 && g[2].abs() <= 1e-12);
    }

    #[test]
    fn lower_exponential_radial_histogram_adds_shell_correction() {
        // Default histogram input on a radial geometry adds the +2kt/ξ shell
        // correction to dU/dξ — the §6.2 Jacobian is load-bearing at lowering.
        let r = lower_one_profile("profile exponential radial 0. 0. 0. lambda 2.5");
        let (lambda, xi): (F, F) = (2.5, 5.0);
        let want_mag = PROFILE_KT / lambda + 2.0 * PROFILE_KT / xi; // 0.4 + 0.4
        let x = [5.0, 0.0, 0.0];
        let mut g = [0.0; 3];
        let _ = r.fg(&x, 1.0, 0.0, &mut g);
        assert!(
            (g[0] - want_mag).abs() <= 1e-9,
            "radial exp+shell force gx: got {}, want {want_mag}",
            g[0]
        );
    }

    #[test]
    fn lower_tabulated_density_radial_is_finite() {
        // A lowered tabulated density evaluates to finite f / fg everywhere on
        // its grid (the §6.4 floor keeps it bounded).
        let r = lower_one_profile(
            "profile tabulated radial 0. 0. 0. density 1. 1.0 3. 0.6 5. 0.2 7. 0.05",
        );
        for &rr in &[2.0, 4.0, 6.0] {
            let x = [rr, 0.0, 0.0];
            let mut g = [0.0; 3];
            let u = r.fg(&x, 1.0, 0.0, &mut g);
            assert!(u.is_finite(), "tabulated U not finite at r={rr}");
            assert!(g.iter().all(|v| v.is_finite()), "tabulated grad not finite");
        }
    }
}
