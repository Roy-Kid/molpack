//! Contract tests for the packing-stage seam (`src/stage.rs`).
//!
//! Everything here runs on **fake stages**. The seam's own behaviour is
//! object safety, the two declaration methods and their composition along a
//! chain, the constructor paths of its three `#[non_exhaustive]` structs, the
//! two default handler hooks, and the re-entrancy contract on
//! [`Stage::run`](molpack::Stage::run) — none of which needs a real algorithm,
//! and all of which a real algorithm would only obscure. What each concrete
//! implementor declares belongs to that implementor's owner
//! (`tests/gencan.rs`, `tests/grow.rs`), so this file boots no engine entry
//! and names no production stage (acceptance ac-008).
//!
//! Fixture: `PackState::new(PackContext::new(0, 0, 0), 0)` — the degenerate
//! context of `tests/geometry_cache.rs` / `src/context/pack_state/tests.rs`,
//! which is all a stage that does no geometry can legitimately need. No RNG,
//! no clock, no filesystem, no network.
//!
//! Single-file gate:
//!
//! ```text
//! cargo test -p molcrafts-molpack --test stage
//! ```

use molpack::handler::StageInfo;
use molpack::{
    Budget, Guarantees, Handler, PackContext, PackState, Placed, Requires, Stage, StageOutcome,
    StepInfo, Target,
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
    softened_per_run: usize,
    /// How many times `run` has been called on this instance.
    runs: usize,
    /// Stand-in for construction-time configuration; `run` must not drain it.
    config: Vec<u32>,
}

impl FakeStage {
    /// A fake declaring `requires` → `guarantees` and reporting
    /// `softened_per_run` relaxations on every run, with a three-element
    /// configuration it is expected to keep.
    fn new(
        name: &'static str,
        requires: Placed,
        guarantees: Placed,
        softened_per_run: usize,
    ) -> Self {
        Self {
            name,
            requires,
            guarantees,
            softened_per_run,
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
    ) -> StageOutcome {
        self.runs += 1;
        // Read the configuration; deliberately do NOT move or drain it — a
        // stage's own configuration must survive every run it is given.
        let _configured: u32 = self.config.iter().sum();
        state.set_placed(self.guarantees);
        StageOutcome::new(true, self.softened_per_run)
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
        outcome.softened, 7,
        "StageOutcome::new takes softened second"
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
/// body is the assertion, exactly as for `NullHandler`'s other defaults.
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

    let first = stage.run(&mut state, &targets, &budget, &mut handlers);
    let config_after_first = stage.config.len();
    assert_eq!(
        state.placed(),
        Placed::All,
        "the first run must move the marker"
    );

    let second = stage.run(&mut state, &targets, &budget, &mut handlers);

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

// ── 8. regression scenario (hard-coded golden) ─────────────────────────────

/// The `Placed` marker after wrapping, after the first stage, after the
/// second. Hard-coded golden, not a computed expectation.
///
/// Provenance: the two fakes' own declarations, written out by hand in this
/// file — `alpha` declares `Placed::None → Placed::All`, `beta` declares
/// `Placed::All → Placed::All`. No tool, no third-party oracle, no real
/// algorithm produced these values; they are the seam's contract spelled as
/// literals so that a change of chaining semantics has to edit them
/// deliberately.
const GOLDEN_PLACED_SEQUENCE: [Placed; 3] = [Placed::None, Placed::All, Placed::All];
/// Hard-coded golden: `alpha` softens twice per run, `beta` three times, and
/// a chain's relaxation count is their sum. `2 + 3 == 5`, written out.
const GOLDEN_SOFTENED_TOTAL: usize = 5;
/// Hard-coded golden: both fakes converge, so a two-stage chain reports two
/// `true`s in order.
const GOLDEN_CONVERGED: [bool; 2] = [true, true];

/// Regression scenario: two fake stages run in sequence on one `PackState`
/// reproduce the hard-coded marker sequence, relaxation sum and convergence
/// flags above (acceptance ac-010).
#[test]
fn stage_regression_fake_chain_outcome_golden() {
    let mut state = empty_state();
    let targets: Vec<Target> = Vec::new();
    let budget = tiny_budget();
    let mut handlers: Vec<Box<dyn Handler>> = Vec::new();

    let mut chain: Vec<Box<dyn Stage>> = vec![
        Box::new(FakeStage::new("alpha", Placed::None, Placed::All, 2)),
        Box::new(FakeStage::new("beta", Placed::All, Placed::All, 3)),
    ];

    let mut observed_placed = vec![state.placed()];
    let mut observed_converged = Vec::with_capacity(chain.len());
    let mut softened_total = 0usize;
    for stage in chain.iter_mut() {
        let outcome = stage.run(&mut state, &targets, &budget, &mut handlers);
        softened_total += outcome.softened;
        observed_converged.push(outcome.converged);
        observed_placed.push(state.placed());
    }

    assert_eq!(
        observed_placed.as_slice(),
        GOLDEN_PLACED_SEQUENCE.as_slice(),
        "the marker transition of a two-stage chain drifted from the golden"
    );
    assert_eq!(
        softened_total, GOLDEN_SOFTENED_TOTAL,
        "a chain's relaxation count is the sum over its stages"
    );
    assert_eq!(
        observed_converged.as_slice(),
        GOLDEN_CONVERGED.as_slice(),
        "both fakes converge, in order"
    );
}
