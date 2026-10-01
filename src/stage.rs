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
//! the pair [`State`](crate::State) reports — are **authoritative
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

pub use crate::outcome::StageOutcome;
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

#[cfg(test)]
mod tests {
    //! Contract tests for the packing-stage seam (`src/stage.rs`).
    //!
    //! Everything here runs on **fake stages**. The seam's own behaviour is
    //! object safety, the two declaration methods and their composition along a
    //! chain, the constructor paths of its three `#[non_exhaustive]` structs, the
    //! two default handler hooks, and the re-entrancy contract on
    //! [`Stage::run`](crate::Stage::run) — none of which needs a real algorithm,
    //! and all of which a real algorithm would only obscure. What each concrete
    //! implementor declares belongs to that implementor's owner
    //! (`gencan::tests`, `grow::tests`), so this file boots no engine entry
    //! and names no production stage (acceptance ac-008).
    //!
    //! Fixture: `PackState::new(PackContext::new(0, 0, 0), 0)` — the degenerate
    //! context of `pack_context.rs::geometry_cache_tests` / `src/context/pack_state/tests.rs`,
    //! which is all a stage that does no geometry can legitimately need. No RNG,
    //! no clock, no filesystem, no network.
    //!
    //! Single-file gate:
    //!
    //! ```text
    //! cargo test -p molcrafts-molpack --lib
    //! ```

    use crate::handler::StageInfo;
    use crate::{
        Budget, Guarantees, Handler, PackContext, PackError, PackState, Placed, Requires, Stage,
        StageOutcome, StepInfo, Target,
    };

    // ── the fake stage ─────────────────────────────────────────────────────────

    /// A stage that moves no atoms.
    ///
    /// It advertises whatever `requires` / `guarantees` markers it is built with,
    /// counts its own `run` calls, and carries a `config` vector standing in for
    /// the construction-time configuration a real stage holds (GENCAN's optimizer
    /// bindings, growth's per-species trees). `run` reads that configuration and
    /// never takes it, which is the re-entrancy contract in miniature.
    struct FakeStage {
        name: &'static str,
        requires: Placed,
        guarantees: Placed,
        degraded_per_run: usize,
        /// How many times `run` has been called on this instance.
        runs: usize,
        /// Stand-in for construction-time configuration; `run` must not drain it.
        config: Vec<u32>,
    }

    impl FakeStage {
        /// A fake declaring `requires` → `guarantees` and reporting
        /// `degraded_per_run` relaxations on every run, with a three-element
        /// configuration it is expected to keep.
        fn new(
            name: &'static str,
            requires: Placed,
            guarantees: Placed,
            degraded_per_run: usize,
        ) -> Self {
            Self {
                name,
                requires,
                guarantees,
                degraded_per_run,
                runs: 0,
                config: vec![11, 22, 33],
            }
        }
    }

    impl Stage for FakeStage {
        fn name(&self) -> &'static str {
            self.name
        }

        fn requires(&self) -> Requires {
            Requires::new(self.requires)
        }

        fn guarantees(&self) -> Guarantees {
            Guarantees::new(self.guarantees)
        }

        fn run(
            &mut self,
            state: &mut PackState,
            _targets: &[Target],
            _budget: &Budget,
            _handlers: &mut [Box<dyn Handler>],
        ) -> Result<StageOutcome, PackError> {
            self.runs += 1;
            // Read the configuration; deliberately do NOT move or drain it — a
            // stage's own configuration must survive every run it is given.
            let _configured: u32 = self.config.iter().sum();
            state.set_placed(self.guarantees);
            Ok(StageOutcome::new(true, self.degraded_per_run))
        }
    }

    /// A handler that implements only the one required method, so every stage
    /// hook it answers is the trait's provided default.
    struct NoopHandler;

    impl Handler for NoopHandler {
        fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {}
    }

    /// The degenerate run state: an empty context and no free molecules.
    fn empty_state() -> PackState {
        PackState::new(PackContext::new(0, 0, 0), 0)
    }

    /// The budget every fake here is run with; no fake reads it.
    fn tiny_budget() -> Budget {
        Budget::new(1, 1e-2)
    }

    // ── 1. object safety ───────────────────────────────────────────────────────

    /// `Stage` must be usable as `Box<dyn Stage>`: a pipeline holds a
    /// heterogeneous, run-time-chosen sequence of stages, so a seam that is not
    /// object safe is not a seam at all.
    #[test]
    fn stage_is_object_safe_behind_a_box() {
        let stages: Vec<Box<dyn Stage>> = vec![
            Box::new(FakeStage::new("alpha", Placed::None, Placed::All, 0)),
            Box::new(FakeStage::new("beta", Placed::All, Placed::All, 0)),
        ];
        let expected = [
            ("alpha", Placed::None, Placed::All),
            ("beta", Placed::All, Placed::All),
        ];

        assert_eq!(stages.len(), expected.len());
        for (stage, &(name, requires, guarantees)) in stages.iter().zip(expected.iter()) {
            assert_eq!(stage.name(), name, "name must be readable through the box");
            assert_eq!(
                stage.requires().placed,
                requires,
                "requires must be readable through the box"
            );
            assert_eq!(
                stage.guarantees().placed,
                guarantees,
                "guarantees must be readable through the box"
            );
        }
    }

    // ── 2. the two declaration methods ─────────────────────────────────────────

    /// A stage reports exactly the markers it was built with: `requires` is the
    /// entry precondition, `guarantees` the exit promise, and neither is derived
    /// from the other.
    #[test]
    fn requires_and_guarantees_report_the_declared_markers() {
        let stage = FakeStage::new("alpha", Placed::None, Placed::All, 0);

        assert_eq!(
            stage.requires().placed,
            Placed::None,
            "a stage that starts from nothing requires Placed::None"
        );
        assert_eq!(
            stage.guarantees().placed,
            Placed::All,
            "a stage that places every free molecule guarantees Placed::All"
        );
    }

    // ── 3. composition along a chain ───────────────────────────────────────────

    /// Two stages compose when the first's exit promise is the second's entry
    /// precondition. The check the pipeline will make is `Placed` equality on
    /// those two markers, and nothing else.
    #[test]
    fn a_chain_links_when_guarantees_meet_the_next_requires() {
        let first = FakeStage::new("alpha", Placed::None, Placed::All, 0);
        let second = FakeStage::new("beta", Placed::All, Placed::All, 0);
        let unlinkable = FakeStage::new("gamma", Placed::None, Placed::All, 0);

        assert_eq!(
            first.guarantees().placed,
            second.requires().placed,
            "alpha guarantees Placed::All, which is exactly what beta requires"
        );
        // Guards the assertion above against being vacuously true: `Placed`
        // equality has to be able to say "no" as well as "yes".
        assert_ne!(
            first.guarantees().placed,
            unlinkable.requires().placed,
            "a successor requiring Placed::None must NOT match a Placed::All promise"
        );
    }

    // ── 4. the constructor paths ───────────────────────────────────────────────

    /// `Requires`, `Guarantees` and `StageOutcome` are `#[non_exhaustive]`, so
    /// each has exactly one way in: its constructor.
    ///
    /// There is no compile-fail harness behind this. An integration test is a
    /// *different* crate from `molcrafts-molpack`, so the struct literals
    /// `Requires { placed: Placed::All }`, `Guarantees { .. }` and
    /// `StageOutcome { .. }` are rejected by the compiler here as a matter of
    /// language rule, not of crate policy — this file could not express them if
    /// it tried, and a harness asserting that the compiler implements
    /// `#[non_exhaustive]` would test rustc, not molpack. What is worth pinning
    /// is the *shape of the constructors*, which is this crate's decision: three
    /// of them, taking the fields in the documented order.
    #[test]
    fn the_three_constructors_fill_the_non_exhaustive_structs() {
        let requires = Requires::new(Placed::All);
        assert_eq!(requires.placed, Placed::All);

        let guarantees = Guarantees::new(Placed::None);
        assert_eq!(guarantees.placed, Placed::None);

        let outcome = StageOutcome::new(true, 7);
        assert!(outcome.converged, "StageOutcome::new takes converged first");
        assert_eq!(
            outcome.degraded, 7,
            "StageOutcome::new takes degraded second"
        );
    }

    // ── 5. the stage identity a handler sees ───────────────────────────────────

    /// `StageInfo` is a plain `Copy` struct like `PhaseInfo`: constructible by
    /// literal, three public fields, all readable. That literal is what lets a
    /// fake — here and in any downstream handler test — drive the two hooks
    /// without booting a pipeline.
    #[test]
    fn stage_info_is_a_literal_with_three_readable_fields() {
        let info = StageInfo {
            index: 0,
            total: 1,
            name: "alpha",
        };

        assert_eq!(
            info.index, 0,
            "the single stage of a one-stage run is index 0"
        );
        assert_eq!(info.total, 1, "a one-stage run reports total 1");
        assert_eq!(info.name, "alpha", "the name is the stage's own name()");
    }

    // ── 6. the two default hooks ───────────────────────────────────────────────

    /// `on_stage_start` / `on_stage_end` are *provided* methods: a handler that
    /// implements only `on_step` still answers both, and answering them does
    /// nothing. The test passes iff neither call panics — reaching the end of the
    /// body is the assertion, exactly as for the trait's other defaults.
    #[test]
    fn the_two_stage_hooks_default_to_no_ops() {
        // Driven through the same `Box<dyn Handler>` shape `Stage::run` receives,
        // so the hooks are pinned on the trait object and not only on the
        // concrete type.
        let mut handler: Box<dyn Handler> = Box::new(NoopHandler);
        let info = StageInfo {
            index: 0,
            total: 1,
            name: "alpha",
        };
        let outcome = StageOutcome::new(false, 0);
        let sys = PackContext::new(0, 0, 0);

        handler.on_stage_start(&info);
        handler.on_stage_end(&info, &outcome, &sys);
    }

    // ── 7. the re-entrancy contract ────────────────────────────────────────────

    /// A stage may be run more than once on an evolving state, and the second run
    /// must have the same capabilities as the first: it may consume the scratch
    /// space it builds per run, never its own configuration.
    ///
    /// The state genuinely evolves between the two calls — the first run advances
    /// the marker from `Placed::None` to `Placed::All` — so this is the multi-run
    /// shape a pipeline produces, not two runs on a pristine state.
    #[test]
    fn a_stage_keeps_its_configuration_across_runs() {
        let mut state = empty_state();
        let targets: Vec<Target> = Vec::new();
        let budget = tiny_budget();
        let mut handlers: Vec<Box<dyn Handler>> = Vec::new();

        let mut stage = FakeStage::new("alpha", Placed::None, Placed::All, 0);
        let config_at_construction = stage.config.len();

        let first = stage
            .run(&mut state, &targets, &budget, &mut handlers)
            .expect("the fake stage runs");
        let config_after_first = stage.config.len();
        assert_eq!(
            state.placed(),
            Placed::All,
            "the first run must move the marker"
        );

        let second = stage
            .run(&mut state, &targets, &budget, &mut handlers)
            .expect("the fake stage runs a second time");

        assert_eq!(stage.runs, 2, "both calls must reach the stage body");
        assert_eq!(
            config_after_first, config_at_construction,
            "the first run consumed part of the stage's configuration"
        );
        assert_eq!(
            stage.config.len(),
            config_at_construction,
            "the second run must find the same configuration the first did"
        );
        assert!(
            first.converged && second.converged,
            "a run that kept its configuration reports the same outcome twice"
        );
    }
}
