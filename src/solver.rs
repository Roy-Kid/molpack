//! The packing-solver seam.
//!
//! [`PackEngine::run`](crate::PackEngine::run) is five stages; the middle
//! two — initial state and the iteration driver — are *the algorithm*, and
//! everything around them (target lowering, `PackContext` construction,
//! frame assembly) is shared infrastructure. A [`Solver`] is one
//! interchangeable implementation of that middle: it receives a fully built
//! context, drives it to a feasible state, and reports the outcome in the
//! same currency (`fdist` / `frest`) every algorithm is judged by.
//!
//! Two rules keep solvers honest peers rather than nested layers:
//!
//! - A solver never calls the GENCAN-internal entry points — those belong
//!   to the rigid-body path, which is itself just another `Solver`
//!   (`gencan::solver::GencanSolver`). The acceptance contract enforces
//!   this with a grep, so even naming them here is off-limits.
//! - The `fdist` / `frest` a solver reports must come from the shared
//!   objective ([`Constraints`](crate::constraints::Constraints)) evaluated
//!   on the final state — never from the solver's own bookkeeping, which may
//!   disagree with the one ruler all algorithms share.
//!
//! **Rust-only:** this module ([`Solver`], [`PlacementsMut`], [`Budget`],
//! [`SolveOutcome`]) is deliberately not mirrored in the Python wheel —
//! Python picks the algorithm by picking the entry (`GenCanPack` /
//! `CbmcGrow`), and implementing a custom solver is a Rust-level extension
//! point.

use crate::context::PackContext;
use crate::handler::Handler;
use crate::target::Target;
use molrs::types::F;

/// One packing algorithm, selected by picking its engine entry
/// ([`GenCanPack`](crate::GenCanPack), [`CbmcGrow`](crate::CbmcGrow)).
pub trait Solver: Send {
    /// Short identifier for logs and reports.
    fn name(&self) -> &'static str;

    /// Drive `sys` to a feasible state.
    ///
    /// `sys` arrives fully built (radii, restraints, `SimBox` + `CellGrid`).
    /// `targets` are the targets this solver is responsible for — the same
    /// objects the caller handed to `pack`, so chemistry has exactly one
    /// source of truth. The solver writes the per-copy conformers into
    /// `sys.coor` and the placements into `x`, and returns the outcome.
    fn solve(
        &mut self,
        sys: &mut PackContext,
        targets: &[Target],
        x: PlacementsMut<'_>,
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> SolveOutcome;
}

/// Typed view over the flat placement vector.
///
/// The layout — COM block first, Euler block second, three values per
/// molecule each (`src/initial.rs::init_xcart_from_x`) — is the GENCAN
/// path's internal convention. Solvers read and write placements through
/// this view and never index the flat vector directly.
#[derive(Debug)]
pub struct PlacementsMut<'a> {
    x: &'a mut [F],
    nmol: usize,
}

impl<'a> PlacementsMut<'a> {
    /// Wrap a flat placement vector holding `6 * nmol` variables.
    ///
    /// # Panics
    ///
    /// Panics when `x.len() != 6 * nmol` — a mismatched view is a programming
    /// error, not a runtime condition.
    pub fn new(x: &'a mut [F], nmol: usize) -> Self {
        assert_eq!(
            x.len(),
            6 * nmol,
            "placement vector holds 6 variables per molecule (3 COM + 3 Euler)"
        );
        Self { x, nmol }
    }

    /// Number of molecules this view addresses.
    pub fn nmol(&self) -> usize {
        self.nmol
    }

    /// Centre of mass of molecule `i`.
    pub fn com(&self, i: usize) -> [F; 3] {
        let o = 3 * i;
        [self.x[o], self.x[o + 1], self.x[o + 2]]
    }

    /// Set the centre of mass of molecule `i`.
    pub fn set_com(&mut self, i: usize, com: [F; 3]) {
        let o = 3 * i;
        self.x[o..o + 3].copy_from_slice(&com);
    }

    /// Euler angles `(beta, gamma, theta)` of molecule `i`.
    pub fn euler(&self, i: usize) -> [F; 3] {
        let o = 3 * self.nmol + 3 * i;
        [self.x[o], self.x[o + 1], self.x[o + 2]]
    }

    /// Set the Euler angles of molecule `i`.
    pub fn set_euler(&mut self, i: usize, euler: [F; 3]) {
        let o = 3 * self.nmol + 3 * i;
        self.x[o..o + 3].copy_from_slice(&euler);
    }

    /// The backing flat vector, for handing to the shared objective.
    pub fn as_slice(&self) -> &[F] {
        self.x
    }

    /// Mutable access to the backing flat vector. The GENCAN path drives
    /// its phase machinery over the raw layout; growth keeps to the typed
    /// accessors above.
    pub fn as_mut_slice(&mut self) -> &mut [F] {
        self.x
    }
}

/// Iteration budget: the engine lifecycle's `max_loops` and `precision`.
/// Each solver documents how it spends the budget.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Budget {
    /// Outer-iteration allowance. The GENCAN path reads this as its loop
    /// count; growth reads it as an allowance of *passes over a chain*, so its
    /// round loop is capped at `max_loops × (the longest species' n_steps + 1)`
    /// rounds — see [`grow::driver`](crate::grow::driver). Serial growth
    /// scheduling advances only one chain per round
    /// ([`GrowConfig::with_serial`](crate::grow::GrowConfig::with_serial)), so
    /// finishing every chain then needs `max_loops ≥ n_chains`.
    pub max_loops: usize,
    /// Convergence threshold on `fdist` / `frest`.
    pub precision: F,
}

impl Budget {
    /// A budget of `max_loops` outer iterations, converged at `precision`.
    pub fn new(max_loops: usize, precision: F) -> Self {
        Self {
            max_loops,
            precision,
        }
    }
}

/// A solver's verdict, fed into [`PackResult`](crate::PackResult).
///
/// `fdist` / `frest` must be produced by the shared objective on the final
/// state (see the module docs) — the seam's guarantee is that every solver
/// is measured with the same ruler.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct SolveOutcome {
    /// Whether the solver reached its own convergence criterion.
    pub converged: bool,
    /// Maximum inter-molecular distance violation on the final state.
    pub fdist: F,
    /// Maximum restraint violation on the final state.
    pub frest: F,
    /// How many times a constructive guarantee had to be relaxed: growth's
    /// hard-core softening rungs plus each forced placement. `0` on the GENCAN
    /// path.
    pub softened: usize,
}

impl SolveOutcome {
    /// A verdict: the solver's own convergence flag, the shared objective's
    /// final `fdist` / `frest`, and how often a guarantee had to be relaxed.
    pub fn new(converged: bool, fdist: F, frest: F, softened: usize) -> Self {
        Self {
            converged,
            fdist,
            frest,
            softened,
        }
    }
}
