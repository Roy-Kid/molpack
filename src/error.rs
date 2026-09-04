use molrs::types::F;
use std::fmt;

#[derive(Debug, Clone)]
pub enum PackError {
    /// Molecules could not satisfy constraints even without distance tolerances.
    ConstraintsFailed(String),
    /// Maximum iterations reached without convergence.
    MaxIterations,
    /// No molecules were provided.
    NoTargets,
    /// A short radius is not shorter than the atom's packing radius, so the
    /// second penalty could never be the tighter one. Packmol rejects the same
    /// input (`app/packmol.f90` lines 519-528).
    ShortRadiusNotShorter {
        target: usize,
        atom: usize,
        short_radius: F,
        radius: F,
    },
    /// A molecule has no atoms.
    EmptyMolecule(usize),
    /// A restraint declared a periodic box whose `max - min` is
    /// non-positive on at least one axis.
    InvalidPBCBox { min: [F; 3], max: [F; 3] },
    /// A declared packing cell is unusable, or contradicts a periodic box.
    InvalidCell { detail: String },
    /// A half-space restraint was declared across a periodic lattice direction,
    /// where it has no well-defined meaning.
    PlaneAcrossPeriodicAxis { axis: usize, normal: [F; 3] },
    /// Two or more restraints declared periodic boxes with different
    /// bounds or different per-axis periodicity flags. Only one periodic
    /// box is allowed per packing run.
    ConflictingPeriodicBoxes {
        first: ([F; 3], [F; 3], [bool; 3]),
        second: ([F; 3], [F; 3], [bool; 3]),
    },
    /// A target was handed to a growth entry (`CbmcGrow`) but cannot be
    /// grown. Growth consumes the template's bond graph; molpack neither
    /// guesses missing chemistry nor silently falls back to rigid-body
    /// packing.
    Grow {
        target: usize,
        source: crate::grow::GrowError,
    },
    /// `with_density` combined with an explicit box or cell — the two are
    /// competing definitions of the same volume.
    DensityConflictsWithBox,
    /// A seeded run's free targets do not match the seed's placement shape
    /// (`GenCanPack::with_restart` — total free atoms expected vs carried).
    SeedMismatch { expected: usize, got: usize },
    /// `with_density` needs every target's mass, and this target's elements
    /// cannot provide one (nor did `Target::with_mass`). Named error, not a
    /// guess.
    UnknownMass { target: usize },
    /// A stage was chained where its entry precondition cannot hold: nothing
    /// before it leaves the placements the stage declares it needs
    /// (`Stage::requires`). Reported before any stage runs and before any
    /// handler is notified.
    StageOrder {
        /// The offending stage's `Stage::name`.
        stage: &'static str,
        /// The precondition that was not met, rendered (e.g. `"placed: all"`).
        needs: &'static str,
    },
    /// A preset entry carrying a non-default *shared* setting was handed to
    /// `Pipeline::with_stage`. The shared settings are one ruler for the whole
    /// run; two stages each holding one would leave the shared objective with
    /// no single ruler, and picking a winner silently is the debt this error
    /// exists to prevent.
    PresetSettingsInsidePipeline {
        /// The name of the first stage that preset produces.
        stage: &'static str,
        /// The `PackSettings` field the preset set (e.g. `"seed"`).
        knob: &'static str,
    },
    /// A pipeline was run with no stages at all.
    NoStages,
    /// A guarded stage finished with an invariant broken, and the guard's
    /// policy was to fail rather than rerun (`OnViolation::Fail`, or a
    /// `Rerun` budget of zero). molpack does not answer this by switching to
    /// another algorithm — the caller picked the method, so the caller picks
    /// the remedy.
    ///
    /// The payload is rendered, deliberately: the layer arrives as its name
    /// (`Layers::name`) and not as a `Layers`, so this module keeps no
    /// dependency on `invariant` or `stage`.
    InvariantViolated {
        /// `Stage::name` of the guarded stage that left the invariant broken.
        stage: &'static str,
        /// `Invariant::name` of the invariant that was broken.
        invariant: &'static str,
        /// The rung of the repair-cost ladder the defect sits on, rendered
        /// (e.g. `"L3 density"`).
        layer: &'static str,
        /// The atoms the violation named, as indices into the run's own
        /// atoms. May be empty when the invariant knows *that* it broke but
        /// not *where*.
        atoms: Vec<usize>,
    },
}

impl fmt::Display for PackError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PackError::ConstraintsFailed(msg) => {
                write!(f, "Packmol failed to satisfy constraints: {msg}")
            }
            PackError::MaxIterations => {
                write!(f, "Maximum iterations reached without convergence")
            }
            PackError::NoTargets => write!(f, "No targets provided"),
            PackError::ShortRadiusNotShorter {
                target,
                atom,
                short_radius,
                radius,
            } => write!(
                f,
                "target {target} atom {atom}: short radius {short_radius} must be \
                 smaller than the packing radius {radius}"
            ),
            PackError::EmptyMolecule(i) => write!(f, "Target {i} has no atoms"),
            PackError::InvalidCell { detail } => write!(f, "invalid packing cell: {detail}"),
            PackError::PlaneAcrossPeriodicAxis { axis, normal } => write!(
                f,
                "plane restraint with normal {normal:?} crosses periodic lattice \
                 direction {axis}: a half-space has no meaning along a periodic axis, \
                 since translating by that lattice vector moves a point across the \
                 plane. Make axis {axis} non-periodic, or drop the plane."
            ),
            PackError::InvalidPBCBox { min, max } => write!(
                f,
                "Invalid PBC box: min={:?}, max={:?} (all max-min components must be > 0)",
                min, max
            ),
            PackError::ConflictingPeriodicBoxes { first, second } => write!(
                f,
                "Conflicting periodic boxes declared by restraints: {first:?} vs {second:?}. \
                 At most one periodic InsideBoxRestraint is allowed per packing run."
            ),
            PackError::Grow { target, source } => {
                write!(f, "target {target} cannot be grown: {source}")
            }
            PackError::DensityConflictsWithBox => write!(
                f,
                "with_density and an explicit box/cell are mutually exclusive: \
                 declare one definition of the volume, not two"
            ),
            PackError::UnknownMass { target } => write!(
                f,
                "with_density needs target {target}'s mass, but its elements do not \
                 resolve one — set Target::with_mass(amu) for element-less species"
            ),
            PackError::StageOrder { stage, needs } => write!(
                f,
                "stage `{stage}` requires {needs} but nothing before it placed the \
                 molecules; put a placing stage (GenCanPack, CbmcGrow, LatticeGrow) \
                 in front of it"
            ),
            PackError::PresetSettingsInsidePipeline { stage, knob } => write!(
                f,
                "preset `{stage}` carries a non-default `{knob}` inside a pipeline; \
                 set `{knob}` on the Pipeline instead (shared settings are one ruler)"
            ),
            PackError::NoStages => write!(
                f,
                "the pipeline has no stages; add one with `with_stage` or run a \
                 preset directly"
            ),
            PackError::InvariantViolated {
                stage,
                invariant,
                layer,
                atoms,
            } => {
                // Bounded rendering: a dense melt can name hundreds of atoms;
                // the count is exact, the list is a preview.
                const SHOWN_ATOMS: usize = 8;
                let shown = &atoms[..atoms.len().min(SHOWN_ATOMS)];
                let more = atoms.len().saturating_sub(SHOWN_ATOMS);
                write!(
                    f,
                    "stage `{stage}` left invariant `{invariant}` ({layer}) broken \
                     on {} atom(s) {shown:?}{}: raise `max` on `OnViolation::Rerun`, \
                     relax the invariant's tolerance, or put a different stage in \
                     front — molpack never switches algorithm for you",
                    atoms.len(),
                    if more > 0 {
                        format!(" and {more} more")
                    } else {
                        String::new()
                    }
                )
            }
            PackError::SeedMismatch { expected, got } => write!(
                f,
                "seeded run: the free targets declare {expected} atoms but the seed \
                 carries {got} — a seeded GenCanPack must receive the same free \
                 targets the seed result was packed from"
            ),
        }
    }
}

impl std::error::Error for PackError {}
