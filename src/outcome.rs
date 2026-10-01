//! What a stage reports about its own successful run.
//!
//! A leaf so the stage seam and the handler can both name it. The run's
//! verdict is not here: it lives on the state the stage just finished
//! writing.

/// What a stage reports about its own *successful* run — and nothing more.
///
/// Both fields are facts only the stage knows. The run's verdict is *not*
/// here: it lives on the state the stage just finished writing, so that
/// every algorithm is measured with the same ruler instead of grading its
/// own paper.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct StageOutcome {
    /// Whether the stage reached its own convergence criterion.
    pub converged: bool,
    /// How many molecules were placed below the stage's own guarantee — one
    /// count per demotion, not per stage.
    ///
    /// `LatticeGrow`: the walk could not keep the occupancy guard (no two
    /// non-bonded atoms closer than the 2nd lattice neighbour), so the chain
    /// went in with site self-avoidance only, or as a forced zigzag. `CbmcGrow`
    /// counts each hard-core softening rung the same way. `0` on the GENCAN
    /// path, which promises nothing constructively and so cannot demote.
    ///
    /// A non-zero count is the honest reading of a crowded box: those
    /// molecules carry the close contacts `fdist` is reporting, and
    /// `converged` is false while it stands.
    pub degraded: usize,
}

impl StageOutcome {
    /// An outcome: the stage's own convergence flag and how often it had to
    /// relax a guarantee.
    pub fn new(converged: bool, degraded: usize) -> Self {
        Self {
            converged,
            degraded,
        }
    }
}
