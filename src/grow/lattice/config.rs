//! Leaf configuration for the lattice growth entry.
//!
//! Same leaf discipline as [`GrowConfig`](crate::grow::config::GrowConfig):
//! only the sibling `prior` module and `molrs` types may be imported here.

use crate::grow::prior::TorsionPrior;

/// Configuration of the diamond-lattice growth solver
/// (lattice-growth-phase spec). The torsion prior is mandatory — it decides
/// the trans/gauche± weights the on-lattice walk grows with — so it is the
/// one constructor argument.
#[derive(Debug, Clone)]
pub struct LatticeConfig {
    pub(crate) torsion_prior: TorsionPrior,
    pub(crate) occupancy_guard: bool,
    pub(crate) max_backtrack: usize,
    pub(crate) max_reseed: usize,
}

impl LatticeConfig {
    pub fn new(torsion_prior: TorsionPrior) -> Self {
        Self {
            torsion_prior,
            occupancy_guard: true,
            max_backtrack: 20_000,
            max_reseed: 200,
        }
    }

    /// Nearest-neighbour site exclusion (on by default): a site may be taken
    /// only when none of its 4 lattice neighbours holds a non-bonded atom,
    /// which keeps every non-bonded pair at ≥ the 2nd-neighbour distance.
    pub fn with_occupancy_guard(mut self, on: bool) -> Self {
        self.occupancy_guard = on;
        self
    }

    /// Recoil budget per chain attempt before the walk reseeds elsewhere.
    pub fn with_max_backtrack(mut self, n: usize) -> Self {
        self.max_backtrack = n.max(1);
        self
    }

    /// Reseed attempts per chain before the solver gives up on the guard
    /// (the escape is counted in `degraded` and the run reports honestly).
    pub fn with_max_reseed(mut self, n: usize) -> Self {
        self.max_reseed = n.max(1);
        self
    }
}
