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
    /// (`GenCanPack::seeded_from` — total free atoms expected vs carried).
    SeedMismatch { expected: usize, got: usize },
    /// `with_density` needs every target's mass, and this target's elements
    /// cannot provide one (nor did `Target::with_mass`). Named error, not a
    /// guess.
    UnknownMass { target: usize },
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
