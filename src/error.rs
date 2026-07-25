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
        }
    }
}

impl std::error::Error for PackError {}
