//! The packing-stage seam.
//!
//! [`PackEngine::run`](crate::PackEngine::run) is five stages; the middle
//! two — initial state and the iteration driver — are *the algorithm*, and
//! everything around them (target lowering, `PackContext` construction,
//! frame assembly) is shared infrastructure. A [`Stage`] is one
//! interchangeable implementation of that middle: it receives the run's
//! [`PackState`], drives it towards a feasible configuration, and declares
//! what it needed on the way in and what it promises on the way out.
//!
//! Every algorithm in this crate is a peer on this seam — a stage never
//! reaches into another stage's driver, and the verdict on a run always
//! comes from the shared objective evaluated on the final state, never from
//! a stage's own bookkeeping.
//!
//! **Rust-only:** this module ([`Stage`], [`Requires`], [`Guarantees`],
//! [`StageOutcome`], [`Budget`]) is deliberately not mirrored in the Python
//! wheel — Python picks the algorithm by picking the entry (`GenCanPack` /
//! `CbmcGrow` / `LatticeGrow`), and implementing a custom stage is a
//! Rust-level extension point.
//!
//! # The repair-cost ladder
//!
//! Structural defects are not equally expensive to repair, and the crate
//! orders them as a six-rung ladder L0–L5. The ladder is a **type**, and it
//! lives with its only reader: [`Layers`](crate::invariant::Layers), next to
//! [`Invariant::layer`](crate::Invariant::layer). Nothing here branches on a
//! rung — this seam's declarations ([`Requires`] / [`Guarantees`]) are about
//! placement shape — so the table and the rung names are documented there,
//! once.
//!
//! **The rule the ladder exists for: a stage is responsible only for the
//! layers it declares.** A stage that promises nothing about chain
//! statistics has not failed when they are poor; a stage that promises no
//! overlaps has failed when overlaps remain. A caller composes a run by
//! stacking stages until every rung it cares about is owned by someone, and
//! guards the ones that matter with
//! [`Pipeline::with_guarded`](crate::Pipeline::with_guarded).
//!
//! # Where the verdict lives
//!
//! The two violation maxima the shared objective produces — the largest
//! inter-molecular contact violation and the largest restraint violation,
//! the pair [`PackResult`](crate::PackResult) reports — are **authoritative
//! on the state after [`Stage::run`] returns**, where the context owns them
//! as its own fields. [`StageOutcome`] carries no verdict: a stage reports
//! only what it alone knows (whether it hit its own convergence criterion,
//! and how many times it had to relax a constructive guarantee). A handler
//! that wants the numbers reads them off the context in
//! [`Handler::on_stage_end`], which is handed the state precisely so that no
//! stage can self-report a verdict the shared ruler would disagree with.
//!
//! # What this seam deliberately does not have
//!
//! * **No `validate` hook.** Not one implementor in this crate would
//!   override it: the rigid-body path validates its targets from its entry,
//!   and both growth paths validate their cell from theirs. A pre-flight
//!   hook nobody implements is a step a caller can forget plus a concept
//!   nobody pays for. The seam is exactly four methods.
//! * **No layer type of its own.** The ladder is
//!   [`Layers`](crate::invariant::Layers), owned by the module that reads it.

use crate::context::{PackState, Placed};
use crate::error::PackError;
use crate::handler::Handler;
use crate::target::Target;
use molrs::types::F;

/// One packing algorithm, selected by picking its engine entry
/// ([`GenCanPack`](crate::GenCanPack), [`CbmcGrow`](crate::CbmcGrow),
/// [`LatticeGrow`](crate::LatticeGrow)).
pub trait Stage: Send {
    /// Short identifier for logs and reports. The same string a
    /// [`StageInfo`](crate::handler::StageInfo) carries to handlers.
    fn name(&self) -> &'static str;

    /// What the state must already hold for this stage to run.
    fn requires(&self) -> Requires;

    /// What the state is promised to hold once this stage returns.
    fn guarantees(&self) -> Guarantees;

    /// Drive `state` towards a feasible configuration.
    ///
    /// The state arrives fully built (radii, restraints, `SimBox` +
    /// `CellGrid`). `targets` are the targets this stage is responsible for
    /// — the same objects the caller handed to
    /// [`PackEngine::run`](crate::PackEngine::run), so chemistry has exactly
    /// one source of truth. A stage takes the context and the rigid
    /// placement vector apart with
    /// [`PackState::rigid_split_mut`], writes the per-copy conformers into
    /// the context's `coor` and the placements into the
    /// [`RigidView`](crate::RigidView), and returns its outcome.
    ///
    /// Those two together are what the run's output is made of: once `run`
    /// returns, the lifecycle rebuilds the lab-frame coordinates from them
    /// with [`RigidView::write_xcart`](crate::RigidView::write_xcart) before
    /// assembling the frame. A stage that works in lab-frame coordinates
    /// directly (both growth drivers do) must therefore capture them back
    /// with
    /// [`RigidView::capture_from_xcart`](crate::RigidView::capture_from_xcart)
    /// before returning; anything left only in the context's `xcart` is
    /// overwritten.
    ///
    /// # Re-entrancy contract
    ///
    /// A `Stage` may be run more than once on an evolving state
    /// (multi-stage pipelines, `Repeat` / `Guarded`). Implementors must not
    /// consume their own configuration: the second `run` must have every
    /// capability the first had. The only thing a run may consume is the
    /// scratch workspace it creates itself.
    ///
    /// The contract is not decorative. A stage that moves its own bound
    /// optimizers, trees or priors out of `self` on the first call keeps
    /// running afterwards — it just runs *degraded*, with no error and no
    /// name for what it lost. Borrow the configuration, do not take it.
    ///
    /// # Failure
    ///
    /// A stage that cannot do its job **fails by a named
    /// [`PackError`]** — it never disguises failure as `converged = false`,
    /// which is the honest report of "I ran and did not reach my criterion"
    /// and nothing else. The lifecycle propagates the error out of
    /// [`PackEngine::run`](crate::PackEngine::run) on the same path as its
    /// own validation errors: the failing stage gets no `on_stage_end`, the
    /// run gets no `on_finish`, and no half-built result is assembled.
    fn run(
        &mut self,
        state: &mut PackState,
        targets: &[Target],
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError>;
}

/// A stage's entry precondition: the placement shape the state must already
/// have.
///
/// One field, because one precondition is all this crate has a reader for.
/// `#[non_exhaustive]`, so [`Requires::new`] is the way in and a second
/// precondition can be added without breaking callers.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Requires {
    /// The placement shape the stage needs on entry.
    pub placed: Placed,
}

impl Requires {
    /// A precondition of `placed`.
    pub fn new(placed: Placed) -> Self {
        Self { placed }
    }
}

/// A stage's exit promise: the placement shape the state has once the stage
/// returns.
///
/// The mirror of [`Requires`], and the marker a caller advances the state
/// with after a run: a chain links when one stage's `Guarantees` meet the
/// next stage's `Requires`.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Guarantees {
    /// The placement shape the stage leaves behind.
    pub placed: Placed,
}

impl Guarantees {
    /// A promise of `placed`.
    pub fn new(placed: Placed) -> Self {
        Self { placed }
    }
}

/// Iteration budget: the engine lifecycle's `max_loops` and `precision`.
/// Each stage documents how it spends the budget.
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
    /// Convergence threshold on the shared objective's violation maxima.
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

/// What a stage reports about its own *successful* run — and nothing more.
///
/// Both fields are facts only the stage knows. The run's verdict is *not*
/// here: it lives on the state the stage just finished writing (see the
/// module docs), so that every algorithm is measured with the same ruler
/// instead of grading its own paper.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct StageOutcome {
    /// Whether the stage reached its own convergence criterion.
    pub converged: bool,
    /// How many times a constructive guarantee had to be relaxed: growth's
    /// hard-core softening rungs plus each forced placement. `0` on the
    /// GENCAN path, which relaxes nothing.
    pub softened: usize,
}

impl StageOutcome {
    /// An outcome: the stage's own convergence flag and how often it had to
    /// relax a guarantee.
    pub fn new(converged: bool, softened: usize) -> Self {
        Self {
            converged,
            softened,
        }
    }
}
