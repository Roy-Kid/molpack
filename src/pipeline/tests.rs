//! Contract tests for the multi-stage lifecycle body (`src/pipeline/`).
//!
//! [`Pipeline`](crate::Pipeline) is the ONE place a molpack run's five
//! phases live — validate, build state, check the stage chain, run each
//! stage, assemble — so the lifecycle's own contract is owned here (law
//! § 11), not spread across the preset entries. What each preset declares
//! stays with that preset; what the seam itself promises stays in
//! [`crate::stage`]'s tests.
//!
//! What this file pins:
//!
//! 1. **Nothing is dropped silently.** A preset's callbacks are adopted by
//!    the pipeline; a preset carrying non-default *shared* settings into
//!    `with_stage` is refused BY NAME; an unmet stage precondition is named
//!    before a single callback fires; an empty pipeline, and a run
//!    with no targets, are named errors rather than no-op runs.
//! 2. **One verdict, one bracket.** `on_start` / `on_finish` bracket the
//!    whole run, `on_stage_start` / `on_stage_end` bracket each stage,
//!    `StepReport.stage` is monotone with `total` = the stage count, and
//!    `degraded` sums across stages.
//! 3. **A combinator is a stage.** `Repeat` runs its body `n` times;
//!    `Guarded` reruns the same stage or fails by name and never switches
//!    algorithm; both report one stage identity however often the body
//!    runs, adopt the body's callbacks and refuse its shared settings.
//!
//! Everything is deterministic by construction — fixed seeds, no wall clock, no
//! filesystem, no network, no third-party oracle.

use crate::EngineSetup;
use crate::callback::{PhaseProgress, StageProgress};
use crate::grow::TorsionPrior;
use crate::test_fixtures::{chain_frame, inside_box};
use crate::{
    Budget, Callback, CbmcGrow, GencanPack, Guarantees, Invariant, Layers, OnViolation,
    PackContext, PackEngine, PackError, PackSettings, PackState, Pipeline, Placed, Requires,
    RestraintsSatisfied, Stage, StageFactory, StageOutcome, State, StepReport, Target, Until,
    Violation,
};
use molrs::op::F;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

/// Rigid water template (positions + packing radii), copied from
/// `pipeline::tests`.
fn water() -> (Vec<[F; 3]>, Vec<F>) {
    (
        vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        vec![1.52, 1.20, 1.20],
    )
}

// ── the four fixtures, each with its shared settings in ONE place ─────────

/// The rigid fixture with a DECLARED periodic box: 30 waters in a 10 Å cell
/// at a 3.5 Å contact tolerance — deliberately over-tight, so one outer loop
/// leaves `converged == false` and `fdist > 0`. That is what makes "the
/// second GENCAN stage continued rather than restarted" observable.
const DENSE_BOX: F = 10.0;
const DENSE_SEED: u64 = 42;
const DENSE_TOL: F = 3.5;
/// One outer loop for the WHOLE pipeline (v1: one budget for every stage).
const DENSE_LOOPS: usize = 1;

fn dense_targets() -> Vec<Target> {
    let (coords, radii) = water();
    vec![Target::from_coords(&coords, &radii, 30).with_name("water")]
}

/// The shared settings every `dense_targets` spelling must agree on.
fn dense_settings<E: PackEngine>(engine: E) -> E {
    engine
        .with_seed(DENSE_SEED)
        .with_tolerance(DENSE_TOL)
        .with_periodic_box([0.0; 3], [DENSE_BOX; 3], [true; 3])
}

/// The rigid fixture with NO box and NO cell declaration: 60 waters held by
/// an ordinary (non-periodic) box restraint. The packing volume is therefore
/// the fall-back box `initial()` synthesizes from `sidemax`, which is the
/// case ac-006 singles out — the GENCAN stage's preamble must install
/// nothing and must not panic.
const FREE_SEED: u64 = 42;
const FREE_TOL: F = 2.0;
const FREE_LOOPS: usize = 20;

fn boxfree_targets() -> Vec<Target> {
    let (coords, radii) = water();
    vec![
        Target::from_coords(&coords, &radii, 60)
            .with_name("water")
            .with_restraint(inside_box([0.0; 3], [14.0; 3])),
    ]
}

fn boxfree_settings<E: PackEngine>(engine: E) -> E {
    engine.with_seed(FREE_SEED).with_tolerance(FREE_TOL)
}

/// The growth fixture, mirroring `grow::tests::seeded_run_contract` so the
/// `CbmcGrow` → `GencanPack::with_restart` comparison is known to be
/// reachable: two 5-bead chains in a generous 20 Å periodic box.
const CHAIN_BOX: F = 20.0;
const CHAIN_SEED: u64 = 9;
const CHAIN_TOL: F = 1.0;
const CHAIN_LOOPS: usize = 60;

fn chain_targets() -> Vec<Target> {
    vec![Target::new(chain_frame(5, 1.5), 2)]
}

fn chain_settings<E: PackEngine>(engine: E) -> E {
    engine
        .with_seed(CHAIN_SEED)
        .with_tolerance(CHAIN_TOL)
        .with_periodic_box([0.0; 3], [CHAIN_BOX; 3], [true; 3])
}

// ── the bitwise ruler ─────────────────────────────────────────────────────

/// Every coordinate, `fdist` and `frest` identical **to the bit**.
///
/// Not a tolerance: the claim under test is that two spellings of the same
/// run take the same arithmetic path, and a tolerance would hide exactly the
/// divergence (a re-`initial()`, a warm geometry cache, a second `radmax`
/// derivation) these tests exist to catch.
fn assert_bitwise_equal(pipeline: &State, direct: &State, what: &str) {
    let (a, b) = (pipeline.positions(), direct.positions());
    assert_eq!(
        a.len(),
        b.len(),
        "{what}: {} atoms from the pipeline vs {} from the direct spelling",
        a.len(),
        b.len()
    );
    for (i, (pa, pb)) in a.iter().zip(b.iter()).enumerate() {
        for (k, (xa, xb)) in pa.iter().zip(pb.iter()).enumerate() {
            assert_eq!(
                xa.to_bits(),
                xb.to_bits(),
                "{what}: atom {i} component {k} is {xa} in the pipeline and \
                 {xb} in the direct spelling — the two spellings of one run \
                 must agree bit for bit, not approximately"
            );
        }
    }
    assert_eq!(
        pipeline.fdist.to_bits(),
        direct.fdist.to_bits(),
        "{what}: fdist {} vs {} — the verdict is read off the same state by \
         both spellings, so it cannot differ",
        pipeline.fdist,
        direct.fdist
    );
    assert_eq!(
        pipeline.frest.to_bits(),
        direct.frest.to_bits(),
        "{what}: frest {} vs {} — a divergence here means the two spellings \
         did not start from the same (cold) geometry cache",
        pipeline.frest,
        direct.frest
    );
}

// ── observers ─────────────────────────────────────────────────────────────

/// Everything a run told its callbacks, in arrival order.
#[derive(Debug, Default)]
struct Tally {
    /// `on_start` calls (must be exactly one per RUN, not per stage).
    starts: usize,
    /// `on_finish` calls (likewise one per run).
    finishes: usize,
    /// `(index, total, name)` of every `on_stage_start`.
    stage_starts: Vec<(usize, usize, &'static str)>,
    /// `(index, total, name)` of every `on_stage_end`.
    stage_ends: Vec<(usize, usize, &'static str)>,
    /// `StepReport.stage` of every `on_step`.
    steps: Vec<(usize, usize, &'static str)>,
}

type Shared = Arc<Mutex<Tally>>;

fn shared() -> Shared {
    Arc::new(Mutex::new(Tally::default()))
}

/// Records every lifecycle callback into a shared [`Tally`], and — when
/// `stop_after` is finite — asks the run to stop once that many `on_step`
/// events have arrived (the `EarlyStopCallback` shape of
/// `grow::tests::Recorder`).
struct Observer {
    tally: Shared,
    stop_after: usize,
}

impl Callback for Observer {
    fn on_start(&mut self, _ntotat: usize, _ntotmol: usize) {
        self.tally.lock().expect("observer mutex").starts += 1;
    }

    fn on_step(&mut self, step: &StepReport, _sys: &PackContext) {
        self.tally.lock().expect("observer mutex").steps.push((
            step.stage.index,
            step.stage.total,
            step.stage.name,
        ));
    }

    fn on_stage_start(&mut self, stage: &StageProgress) {
        self.tally
            .lock()
            .expect("observer mutex")
            .stage_starts
            .push((stage.index, stage.total, stage.name));
    }

    fn on_stage_end(&mut self, stage: &StageProgress, _outcome: &StageOutcome, _sys: &PackContext) {
        self.tally.lock().expect("observer mutex").stage_ends.push((
            stage.index,
            stage.total,
            stage.name,
        ));
    }

    fn on_finish(&mut self, _sys: &PackContext) {
        self.tally.lock().expect("observer mutex").finishes += 1;
    }

    fn should_stop(&self) -> bool {
        self.tally.lock().expect("observer mutex").steps.len() >= self.stop_after
    }
}

/// An observer that never asks for a stop.
fn observer(tally: &Shared) -> Box<dyn Callback> {
    Box::new(Observer {
        tally: Arc::clone(tally),
        stop_after: usize::MAX,
    })
}

// ── a stage that cannot legally run first ─────────────────────────────────

/// A stage that demands `Placed::All` on entry and panics if it is ever run.
///
/// The panic is the assertion: the stage-order check must reject the chain
/// BEFORE any stage executes, so reaching `run` at all is the failure.
struct NeedsPlacedStage;

impl Stage for NeedsPlacedStage {
    fn name(&self) -> &'static str {
        "needs-placed"
    }

    fn requires(&self) -> Requires {
        Requires::new(Placed::All)
    }

    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    fn run(
        &mut self,
        _state: &mut PackState,
        _targets: &[Target],
        _budget: &Budget,
        _callbacks: &mut [Box<dyn Callback>],
    ) -> Result<StageOutcome, PackError> {
        panic!(
            "the stage-order check let a stage requiring Placed::All run first \
             — the chain must be rejected before any stage executes"
        );
    }
}

/// The minimal [`StageFactory`] a caller can write outside this crate: it
/// carries default shared settings and produces the one fake stage.
struct NeedsPlacedFactory {
    settings: PackSettings,
}

impl NeedsPlacedFactory {
    fn new() -> Self {
        Self {
            settings: PackSettings::default(),
        }
    }
}

impl StageFactory for NeedsPlacedFactory {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    fn stages(&mut self, _setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        Ok(vec![Box::new(NeedsPlacedStage)])
    }
}

// ── 1. one stage in a pipeline IS the preset run ──────────────────────────

// ── 2. cross-algorithm hand-off ≡ with_restart ─────────────────────────────

// ── 3. callbacks are adopted, never dropped ────────────────────────────────

/// A callback attached to a PRESET that is then handed to `with_stage` still
/// receives its callbacks: the pipeline adopts it (ac-005). Dropping it is
/// the silent failure this whole design exists to prevent.
#[test]
fn pipeline_adopts_preset_callbacks() {
    let one = shared();
    let piped = boxfree_settings(
        Pipeline::new().with_stage(GencanPack::new().with_callback(observer(&one))),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a single-stage pipeline with an adopted callback runs");
    assert!(piped.converged);
    {
        let t = one.lock().expect("observer mutex");
        assert!(
            !t.steps.is_empty(),
            "the callback carried in by GencanPack::with_callback saw no \
             on_step events — with_stage must adopt a preset's callbacks, not \
             drop them on the floor"
        );
        assert_eq!(
            t.starts, 1,
            "an adopted callback is bracketed by the run like any other"
        );
        assert_eq!(t.finishes, 1);
    }

    // Two stages: an adopted callback observes the WHOLE run, so both stage
    // indices show up on the events it recorded.
    let both = shared();
    dense_settings(
        Pipeline::new()
            .with_stage(GencanPack::new().with_callback(observer(&both)))
            .with_stage(GencanPack::new()),
    )
    .run(&dense_targets(), DENSE_LOOPS)
    .expect("a two-stage pipeline with an adopted callback runs");

    let t = both.lock().expect("observer mutex");
    let mut seen: Vec<usize> = t.steps.iter().map(|&(index, _, _)| index).collect();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(
        seen,
        vec![0, 1],
        "an adopted callback must see BOTH stages (recorded stage indices \
         {seen:?}) — adoption is for the run, not for the stage that carried \
         it in"
    );
}

// ── 4. the two named rejections and the empty pipeline ────────────────────

/// An unmet `requires()` is reported by name BEFORE any callback is notified
/// and before any stage runs (ac-003).
#[test]
fn pipeline_stage_order_error_fires_before_any_callback() {
    let tally = shared();
    let err = boxfree_settings(
        Pipeline::new()
            .with_callback(observer(&tally))
            .with_stage(NeedsPlacedFactory::new())
            .with_stage(GencanPack::new()),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect_err("a stage requiring Placed::All cannot be first");

    let (stage, needs) = match err {
        PackError::StageOrder { stage, needs } => (stage, needs),
        other => panic!("expected PackError::StageOrder, got {other:?}"),
    };
    assert_eq!(
        stage, "needs-placed",
        "the error must name the offending stage"
    );
    assert!(
        !needs.is_empty(),
        "the error must render the precondition that was not met"
    );
    let msg = format!(
        "{}",
        PackError::StageOrder {
            stage: "needs-placed",
            needs
        }
    );
    assert!(
        msg.contains(stage) && msg.contains(needs),
        "Display must name the stage and the missing precondition, got: {msg}"
    );

    let t = tally.lock().expect("observer mutex");
    assert_eq!(
        t.starts, 0,
        "on_start fired before the chain was checked — the order check runs \
         before ANY callback is notified"
    );
    assert!(
        t.stage_starts.is_empty() && t.stage_ends.is_empty() && t.steps.is_empty(),
        "a rejected chain must produce no stage callbacks at all, got \
         {} stage starts / {} stage ends / {} steps",
        t.stage_starts.len(),
        t.stage_ends.len(),
        t.steps.len()
    );
    assert_eq!(t.finishes, 0, "a rejected chain never finishes a run");
}

/// A preset carrying a non-default SHARED setting into `with_stage` is
/// refused by the name of the knob (ac-004): two stages each holding a ruler
/// would leave the shared objective with no single ruler, and picking a
/// winner silently is exactly the debt this error exists to prevent.
///
/// The error surfaces from `run`, not from `with_stage` — `with_stage`
/// returns `Self`, so it has nowhere to put a `Result`.
#[test]
fn pipeline_rejects_preset_with_non_default_settings_by_name() {
    let err = Pipeline::new()
        .with_stage(GencanPack::new().with_seed(7))
        .run(&boxfree_targets(), 1)
        .expect_err("a preset carrying a seed into with_stage must be refused");

    match err {
        PackError::PresetSettingsInsidePipeline { stage, knob } => {
            assert_eq!(stage, "gencan", "the error names the offending stage");
            assert_eq!(
                knob, "seed",
                "the error names the knob the preset set, so the user knows \
                 what to move onto the pipeline"
            );
            let msg = format!(
                "{}",
                PackError::PresetSettingsInsidePipeline { stage, knob }
            );
            assert!(
                msg.contains(stage) && msg.contains(knob),
                "Display must name stage and knob, got: {msg}"
            );
        }
        other => panic!("expected PackError::PresetSettingsInsidePipeline, got {other:?}"),
    }

    // The same preset run directly (i.e. through `Pipeline::single`) is NOT
    // an error: `single` adopts the settings instead of refusing them.
    let direct = GencanPack::new()
        .with_seed(7)
        .with_tolerance(FREE_TOL)
        .run(&boxfree_targets(), 1)
        .expect("a preset's own run adopts its settings and must not error");
    assert!(direct.fdist.is_finite() && direct.frest.is_finite());

    // …and neither is the explicit `Pipeline::single` spelling.
    let single = Pipeline::single(GencanPack::new().with_seed(7).with_tolerance(FREE_TOL))
        .run(&boxfree_targets(), 1)
        .expect("Pipeline::single adopts the engine's settings");
    assert_bitwise_equal(&single, &direct, "Pipeline::single vs the preset run");
}

/// A pipeline with no stages is a named error, never a run that quietly
/// produces the input back.
#[test]
fn pipeline_empty_is_a_named_error() {
    let err = Pipeline::new()
        .with_seed(FREE_SEED)
        .run(&boxfree_targets(), FREE_LOOPS)
        .expect_err("an empty pipeline has nothing to run");
    assert!(
        matches!(err, PackError::NoStages),
        "expected PackError::NoStages, got {err:?}"
    );
    let msg = format!("{err}");
    assert!(
        msg.contains("stage"),
        "Display must say the pipeline has no stages, got: {msg}"
    );
}

/// No targets is the other named refusal: a run with nothing to place is an
/// error, not an empty success.
#[test]
fn pipeline_without_targets_is_a_named_error() {
    let err = Pipeline::single(GencanPack::new())
        .run(&[], FREE_LOOPS)
        .expect_err("a run with no targets has nothing to place");
    assert!(
        matches!(err, PackError::NoTargets),
        "expected PackError::NoTargets, got {err:?}"
    );
}

// ── 5. one verdict, one bracket per run, one bracket per stage ────────────

/// A two-stage run brackets the RUN once and each STAGE once, reports
/// `stage.total == 2` on every step with a non-decreasing `stage.index`, and
/// sums `degraded` across the stages (ac-009).
#[test]
fn pipeline_two_stages_sum_degraded_and_count_hooks() {
    let targets = chain_targets();

    // The two standalone spellings, same seed and settings, provide the
    // reference sum.
    let grown = chain_settings(CbmcGrow::new(TorsionPrior::Uniform))
        .run(&targets, CHAIN_LOOPS)
        .expect("the growth stage runs standalone");
    let seeded = GencanPack::new()
        .with_restart(&grown)
        .with_seed(CHAIN_SEED)
        .with_tolerance(CHAIN_TOL)
        .run(&targets, CHAIN_LOOPS)
        .expect("the seeded GENCAN stage runs standalone");

    let carried = shared();
    let attached = shared();
    let piped = chain_settings(
        Pipeline::new()
            .with_callback(observer(&attached))
            .with_stage(CbmcGrow::new(TorsionPrior::Uniform).with_callback(observer(&carried)))
            .with_stage(GencanPack::new()),
    )
    .run(&targets, CHAIN_LOOPS)
    .expect("a growth → GENCAN pipeline runs");

    assert_eq!(
        piped.degraded,
        grown.degraded + seeded.degraded,
        "the pipeline's degraded count is the SUM over its stages \
         ({} + {}), not the last stage's own count",
        grown.degraded,
        seeded.degraded
    );

    for (label, tally) in [("attached", &attached), ("adopted", &carried)] {
        let t = tally.lock().expect("observer mutex");
        assert_eq!(
            t.starts, 1,
            "{label}: on_start is once per RUN, not once per stage"
        );
        assert_eq!(t.finishes, 1, "{label}: on_finish is once per RUN");
        assert_eq!(
            t.stage_starts.len(),
            2,
            "{label}: on_stage_start fires once per stage, got {:?}",
            t.stage_starts
        );
        assert_eq!(
            t.stage_ends.len(),
            2,
            "{label}: on_stage_end fires once per stage, got {:?}",
            t.stage_ends
        );
        let indices: Vec<usize> = t.stage_starts.iter().map(|&(i, _, _)| i).collect();
        assert_eq!(
            indices,
            vec![0, 1],
            "{label}: stage indices arrive in chain order"
        );
        let names: Vec<&'static str> = t.stage_starts.iter().map(|&(_, _, n)| n).collect();
        assert_eq!(
            names,
            vec!["growth", "gencan"],
            "{label}: each StageProgress carries the stage's own name()"
        );
        for &(index, total, _) in t.stage_starts.iter().chain(t.stage_ends.iter()) {
            assert_eq!(
                total, 2,
                "{label}: StageProgress.total is the number of stages in the run"
            );
            assert!(index < total, "{label}: stage index {index} out of range");
        }
        for &(index, total, _) in t.steps.iter() {
            assert_eq!(
                total, 2,
                "{label}: every StepReport reports the pipeline's stage count"
            );
            assert!(index < total);
        }
        for pair in t.steps.windows(2) {
            let (prev, cur) = (pair[0].0, pair[1].0);
            assert!(
                cur >= prev,
                "{label}: StepReport.stage.index went {prev} → {cur} — the \
                 stage index is monotone along a linear pipeline"
            );
        }
    }
}

// ── 6. the box-free single-stage run ──────────────────────────────────────

// ── 8. the two combinators (spec 06) ──────────────────────────────────────
//
// `Repeat` and `Guarded` are stages, so everything above still applies to
// them: they are resolved by the same chain check, bracketed by the same
// callbacks, and report ONE stage identity no matter how often their body
// runs. What is new here is the body's run count, the honest `degraded` sum,
// the continuation of a repeated GENCAN pass, and the named refusal a broken
// invariant produces.

/// A stage that moves no atoms and counts its runs through a shared counter.
///
/// The combinators are about *how often* and *under what condition* a stage
/// runs, so the body they wrap only has to be countable — a real algorithm
/// would add arithmetic that hides the count. When asked, the stage also
/// emits one `on_phase_start` per run: that is how a callback learns the body
/// did a unit of work from *inside* a combinator.
struct CountingStage {
    runs: Arc<AtomicUsize>,
    converged: bool,
    degraded: usize,
    signal: bool,
}

impl Stage for CountingStage {
    fn name(&self) -> &'static str {
        "counting"
    }

    fn requires(&self) -> Requires {
        Requires::new(Placed::None)
    }

    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    fn run(
        &mut self,
        _state: &mut PackState,
        _targets: &[Target],
        _budget: &Budget,
        callbacks: &mut [Box<dyn Callback>],
    ) -> Result<StageOutcome, PackError> {
        self.runs.fetch_add(1, Ordering::Relaxed);
        if self.signal {
            for h in callbacks.iter_mut() {
                h.on_phase_start(&PhaseProgress {
                    phase: 0,
                    total_phases: 1,
                    molecule_type: None,
                });
            }
        }
        Ok(StageOutcome::new(self.converged, self.degraded))
    }
}

/// The factory that produces one [`CountingStage`], carrying default shared
/// settings so a pipeline has nothing to refuse.
struct CountingFactory {
    settings: PackSettings,
    runs: Arc<AtomicUsize>,
    converged: bool,
    degraded: usize,
    signal: bool,
}

impl CountingFactory {
    fn new(runs: &Arc<AtomicUsize>) -> Self {
        Self {
            settings: PackSettings::default(),
            runs: Arc::clone(runs),
            converged: false,
            degraded: 0,
            signal: false,
        }
    }

    /// The produced stage reports its own convergence criterion as met.
    fn converging(mut self) -> Self {
        self.converged = true;
        self
    }

    /// The produced stage reports `n` relaxations per run.
    fn with_degraded(mut self, n: usize) -> Self {
        self.degraded = n;
        self
    }

    /// The produced stage emits one `on_phase_start` per run.
    fn signalling(mut self) -> Self {
        self.signal = true;
        self
    }
}

impl StageFactory for CountingFactory {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    fn stages(&mut self, _setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        Ok(vec![Box::new(CountingStage {
            runs: Arc::clone(&self.runs),
            converged: self.converged,
            degraded: self.degraded,
            signal: self.signal,
        })])
    }
}

/// One factory as a combinator body.
fn body(factory: impl StageFactory + 'static) -> Vec<Box<dyn StageFactory>> {
    vec![Box::new(factory)]
}

/// An invariant that is never satisfied — the shape `OnViolation` is defined
/// against. It reads nothing from the state, so the rerun count it produces
/// is the combinator's own arithmetic and not a fixture's.
struct AlwaysViolated;

impl Invariant for AlwaysViolated {
    fn name(&self) -> &'static str {
        "always-violated"
    }

    fn layer(&self) -> Layers {
        Layers::L4_LOCAL_OVERLAPS
    }

    fn check(&self, _state: &PackState) -> Vec<Violation> {
        vec![Violation {
            atoms: vec![0],
            what: "this invariant never holds".to_string(),
        }]
    }
}

/// A callback that asks for a stop as soon as the body has signalled one unit
/// of work. `should_stop` must be honoured from *inside* a combinator, not
/// only between the pipeline's own stages.
struct StopOnFirstSignal {
    signals: Arc<AtomicUsize>,
}

impl Callback for StopOnFirstSignal {
    fn on_step(&mut self, _step: &StepReport, _sys: &PackContext) {}

    fn on_phase_start(&mut self, _phase: &PhaseProgress) {
        self.signals.fetch_add(1, Ordering::Relaxed);
    }

    fn should_stop(&self) -> bool {
        self.signals.load(Ordering::Relaxed) >= 1
    }
}

/// `Until::Passes(n)` runs the body exactly `n` times and sums what each pass
/// relaxed — the same honest accumulation a linear chain performs, so a
/// repeated stage cannot under-report by reusing the last pass's count.
#[test]
fn repeat_passes_runs_body_n_times_and_sums_degraded() {
    const DEGRADED_PER_RUN: usize = 3;
    let runs = Arc::new(AtomicUsize::new(0));

    let result = boxfree_settings(Pipeline::new().with_repeat(
        body(CountingFactory::new(&runs).with_degraded(DEGRADED_PER_RUN)),
        Until::Passes(2),
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a Repeat{Passes(2)} pipeline runs");

    assert_eq!(
        runs.load(Ordering::Relaxed),
        2,
        "Until::Passes(2) must run the body exactly twice"
    );
    assert_eq!(
        result.degraded,
        2 * DEGRADED_PER_RUN,
        "the repeated stage's degraded count is the SUM over its passes \
         (2 × {DEGRADED_PER_RUN}), not one pass's own count"
    );
}

/// `Until::Converged` stops after the first pass whose body converged: the
/// body runs once, not `n` times, and never "one more for luck".
#[test]
fn repeat_until_converged_stops_after_first_converged_pass() {
    let runs = Arc::new(AtomicUsize::new(0));

    boxfree_settings(Pipeline::new().with_repeat(
        body(CountingFactory::new(&runs).converging()),
        Until::Converged,
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a Repeat{Converged} pipeline runs");

    assert_eq!(
        runs.load(Ordering::Relaxed),
        1,
        "the body converged on its first pass, so Until::Converged must not \
         run it again"
    );
}

/// `Until::Passes(0)` contributes NO stage — not a silently clamped one pass.
/// A pipeline holding only that is therefore empty, which is the named
/// [`PackError::NoStages`]; behind a real stage it is simply a no-op.
#[test]
fn repeat_passes_zero_contributes_no_stage() {
    let alone = Arc::new(AtomicUsize::new(0));
    let err = boxfree_settings(
        Pipeline::new().with_repeat(body(CountingFactory::new(&alone)), Until::Passes(0)),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect_err("a pipeline whose only combinator contributes no stage is empty");
    assert!(
        matches!(err, PackError::NoStages),
        "expected PackError::NoStages, got {err:?}"
    );
    assert_eq!(
        alone.load(Ordering::Relaxed),
        0,
        "Passes(0) must not run the body at all"
    );

    let after = Arc::new(AtomicUsize::new(0));
    let result = boxfree_settings(
        Pipeline::new()
            .with_stage(GencanPack::new())
            .with_repeat(body(CountingFactory::new(&after)), Until::Passes(0)),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a pipeline with one real stage runs even if a combinator is empty");
    assert!(result.fdist.is_finite() && result.frest.is_finite());
    assert_eq!(
        after.load(Ordering::Relaxed),
        0,
        "Passes(0) after a real stage still runs the body zero times"
    );
}

/// `OnViolation::Fail` surfaces a NAMED error out of `Pipeline::run`, naming
/// the stage, the invariant, its rung on the repair-cost ladder and the atoms
/// involved — never a quietly unconverged result and never another algorithm
/// (law P8).
///
/// The atom list itself is exercised on a hand-built state in
/// `invariant::tests`: `frest_atom` is movebad-only bookkeeping
/// (`src/movebad.rs`), so which atoms a real run's final state attributes the
/// residual to is the GENCAN schedule's business, not this file's.
#[test]
fn guarded_fail_returns_named_error() {
    // frest ≈ 3.0e-4 on this fixture, so a tolerance of exactly 0.0 is
    // unsatisfiable by construction.
    let invariants: Vec<Box<dyn Invariant>> = vec![Box::new(RestraintsSatisfied::new(0.0))];

    let err = boxfree_settings(Pipeline::new().with_guarded(
        GencanPack::new(),
        invariants,
        OnViolation::Fail,
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect_err("an unsatisfiable guard under OnViolation::Fail must fail by name");

    match &err {
        PackError::InvariantViolated {
            stage,
            invariant,
            layer,
            atoms,
        } => {
            assert_eq!(*stage, "gencan", "the error names the guarded stage");
            assert!(
                !invariant.is_empty(),
                "the error must name the invariant that was broken"
            );
            assert!(
                !layer.is_empty(),
                "the error must render the layer the defect sits on"
            );
            for (k, &icart) in atoms.iter().enumerate() {
                assert!(
                    icart < 180,
                    "reported atom {k} is index {icart}, outside the \
                     fixture's 60 waters × 3 atoms — the payload must index \
                     the run's own atoms"
                );
            }
            let msg = format!("{err}");
            assert!(
                msg.contains(*invariant) && msg.contains(*layer),
                "Display must name the invariant and its layer, got: {msg}"
            );
        }
        other => panic!("expected PackError::InvariantViolated, got {other:?}"),
    }
}

/// `OnViolation::Rerun { max }` reruns the SAME stage at most `max` more
/// times and then reports honestly: `Ok` with `converged == false`, never a
/// dead loop and never a switch to another algorithm. `Rerun { max: 0 }` is
/// `Fail`.
#[test]
fn guarded_rerun_max_two_reruns_twice_then_unconverged() {
    let runs = Arc::new(AtomicUsize::new(0));
    let invariants: Vec<Box<dyn Invariant>> = vec![Box::new(AlwaysViolated)];

    // The body reports `converged == true`; the guard's verdict must win.
    let result = boxfree_settings(Pipeline::new().with_guarded(
        CountingFactory::new(&runs).converging(),
        invariants,
        OnViolation::Rerun { max: 2 },
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("an exhausted Rerun budget is a result, not an error");

    assert_eq!(
        runs.load(Ordering::Relaxed),
        3,
        "Rerun{{max: 2}} is one run plus at most two reruns of the SAME stage"
    );
    assert!(
        !result.converged,
        "a stage that failed its guard on every attempt must not be reported \
         as converged, whatever the stage itself claims"
    );

    let zero_runs = Arc::new(AtomicUsize::new(0));
    let zero_invariants: Vec<Box<dyn Invariant>> = vec![Box::new(AlwaysViolated)];
    let err = boxfree_settings(Pipeline::new().with_guarded(
        CountingFactory::new(&zero_runs).converging(),
        zero_invariants,
        OnViolation::Rerun { max: 0 },
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect_err("Rerun{max: 0} has no budget left, so it is Fail");
    match &err {
        PackError::InvariantViolated {
            stage, invariant, ..
        } => {
            assert_eq!(*stage, "counting");
            assert_eq!(*invariant, "always-violated");
        }
        other => panic!("expected PackError::InvariantViolated, got {other:?}"),
    }
    assert_eq!(
        zero_runs.load(Ordering::Relaxed),
        1,
        "Rerun{{max: 0}} runs the stage once and then fails, exactly like Fail"
    );
}

/// A combinator is ONE stage to the run around it: the inner passes never
/// move `StepReport.stage.index` or `total`, and the names a callback sees are
/// the top-level ones. Nesting a combinator inside another does not change
/// that (ac-006).
#[test]
fn combinators_keep_stage_index_monotone() {
    let flat = shared();
    dense_settings(
        Pipeline::new()
            .with_callback(observer(&flat))
            .with_stage(GencanPack::new())
            .with_repeat(body(GencanPack::new()), Until::Passes(2)),
    )
    .run(&dense_targets(), DENSE_LOOPS)
    .expect("a [gencan, repeat] pipeline runs");

    {
        let t = flat.lock().expect("observer mutex");
        let names: Vec<&'static str> = t.stage_starts.iter().map(|&(_, _, n)| n).collect();
        assert_eq!(
            names,
            vec!["gencan", "repeat"],
            "a combinator reports itself as ONE stage, under its own name"
        );
        for &(index, total, _) in t.steps.iter() {
            assert_eq!(
                total, 2,
                "the run has two TOP-LEVEL stages; an inner pass must not \
                 change the total"
            );
            assert!(index < total, "stage index {index} out of range");
        }
        for pair in t.steps.windows(2) {
            let (prev, cur) = (pair[0].0, pair[1].0);
            assert!(
                cur >= prev,
                "StepReport.stage.index went {prev} → {cur} — a repeated body \
                 must not walk the index backwards"
            );
        }
    }

    // Nested: Guarded(Repeat([gencan])). The invariant list is empty on
    // purpose — this test is about stage identity, and the guard's own
    // verdict is pinned by the two tests above.
    let nested = shared();
    let inner = Pipeline::new().with_repeat(body(GencanPack::new()), Until::Passes(2));
    dense_settings(
        Pipeline::new()
            .with_callback(observer(&nested))
            .with_guarded(inner, Vec::new(), OnViolation::Fail),
    )
    .run(&dense_targets(), DENSE_LOOPS)
    .expect("a Guarded(Repeat([gencan])) pipeline runs");

    let t = nested.lock().expect("observer mutex");
    let names: Vec<&'static str> = t.stage_starts.iter().map(|&(_, _, n)| n).collect();
    assert_eq!(
        names,
        vec!["guarded"],
        "the outer combinator is the only top-level stage"
    );
    for &(index, total, _) in t.steps.iter() {
        assert_eq!(total, 1, "one top-level stage means total == 1");
        assert_eq!(index, 0, "…and index 0, however deep the nesting goes");
    }
}

/// A body factory's callbacks are **adopted** exactly as `with_stage` adopts
/// them, and its non-default shared settings are refused BY NAME exactly as
/// `with_stage` refuses them — a combinator is not a hole in either rule.
#[test]
fn repeat_adopts_body_callbacks_and_refuses_settings() {
    let tally = shared();
    boxfree_settings(Pipeline::new().with_repeat(
        body(GencanPack::new().with_callback(observer(&tally))),
        Until::Passes(2),
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a Repeat whose body carries a callback runs");

    {
        let t = tally.lock().expect("observer mutex");
        assert!(
            !t.steps.is_empty(),
            "the callback carried in by the body preset saw no on_step events \
             — with_repeat must adopt a body's callbacks, not drop them"
        );
        assert_eq!(t.starts, 1, "an adopted callback is bracketed once per RUN");
        assert_eq!(t.finishes, 1);
    }

    let err = Pipeline::new()
        .with_repeat(body(GencanPack::new().with_seed(7)), Until::Passes(2))
        .run(&boxfree_targets(), 1)
        .expect_err("a body preset carrying a seed must be refused");
    match &err {
        PackError::PresetSettingsInsidePipeline { stage, knob } => {
            assert_eq!(
                *knob, "seed",
                "the refusal names the knob so the user knows what to move \
                 onto the pipeline"
            );
            let msg = format!("{err}");
            assert!(
                msg.contains(*stage) && msg.contains(*knob),
                "Display must name stage and knob, got: {msg}"
            );
        }
        other => panic!("expected PackError::PresetSettingsInsidePipeline, got {other:?}"),
    }
}

/// A callback asking to stop from inside a combinator body ends the run there:
/// the remaining passes never start and the verdict is honestly
/// `converged == false`.
#[test]
fn repeat_should_stop_yields_immediately() {
    let runs = Arc::new(AtomicUsize::new(0));
    let signals = Arc::new(AtomicUsize::new(0));

    let result = boxfree_settings(
        Pipeline::new()
            .with_callback(Box::new(StopOnFirstSignal {
                signals: Arc::clone(&signals),
            }))
            .with_repeat(
                body(CountingFactory::new(&runs).converging().signalling()),
                Until::Passes(3),
            ),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a run stopped from inside a combinator still returns Ok");

    assert_eq!(
        runs.load(Ordering::Relaxed),
        1,
        "the body signalled a stop during its first pass, so passes 2 and 3 \
         must never start"
    );
    assert!(
        !result.converged,
        "a run abandoned mid-combinator must never be reported as converged"
    );
}
