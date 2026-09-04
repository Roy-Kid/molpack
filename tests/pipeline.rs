//! Contract tests for the multi-stage lifecycle body
//! (`.claude/specs/stage-pipeline-05-pipeline.md`,
//! `.claude/specs/stage-pipeline-06-combinators.md`, `src/pipeline/`).
//!
//! [`Pipeline`](molpack::Pipeline) is the ONE place a molpack run's five
//! phases live — validate, build state, check the stage chain, run each
//! stage, assemble — so the lifecycle's own contract is owned here (law
//! § 11), not spread across the three preset entries. What each preset
//! declares stays with that preset (`tests/gencan.rs`, `tests/grow.rs`,
//! `tests/packer.rs`); what the seam itself promises stays in
//! `tests/stage.rs`.
//!
//! Three things this file pins that nothing else can:
//!
//! 1. **Two spellings, one answer.** A single-stage pipeline is bitwise the
//!    preset run (ac-006); `[CbmcGrow, GenCanPack]` is bitwise
//!    `GenCanPack::with_restart(cbmc_result)` (ac-007); `[GenCanPack,
//!    GenCanPack]` is bitwise `GenCanPack::with_restart(first_result)`
//!    (ac-008) — the last of which is the observable form of "a second
//!    GENCAN stage continues from the placements it was handed and never
//!    re-runs `initial()`".
//! 2. **Nothing is dropped silently.** A preset's handlers are adopted by
//!    the pipeline (ac-005); a preset carrying non-default *shared*
//!    settings into `with_stage` is refused BY NAME (ac-004); an unmet
//!    stage precondition is named before a single handler callback fires
//!    (ac-003); an empty pipeline is a named error, not a no-op run.
//! 3. **One verdict, honest early stop.** `on_start` / `on_finish` bracket
//!    the whole run, `on_stage_start` / `on_stage_end` bracket each stage,
//!    `StepInfo.stage` is monotone with `total` = the stage count, and a
//!    mid-run `should_stop` leaves `converged == false` with the growth
//!    stage's bonded geometry intact (ac-009).
//! 4. **A combinator is a stage.** `Repeat` runs its body `n` times and its
//!    second pass CONTINUES the first (bitwise `with_restart`, 06 ac-002);
//!    `Guarded` reruns the same stage or fails by name and never switches
//!    algorithm (06 ac-003); both report one stage identity however often
//!    the body runs (06 ac-006), adopt the body's handlers and refuse its
//!    shared settings.
//!
//! Fixtures are copied, never `mod`-shared: the chain / water templates
//! below mirror `tests/grow.rs` and `tests/packer.rs` so that a change to
//! either file cannot silently move this file's answers. Everything is
//! deterministic by construction — fixed seeds, no wall clock, no
//! filesystem, no network, no third-party oracle.
//!
//! Single-file gate:
//!
//! ```text
//! cargo test -p molcrafts-molpack --test pipeline
//! ```

use molpack::grow::TorsionPrior;
use molpack::handler::StageInfo;
use molpack::pipeline::EngineSetup;
use molpack::{
    Budget, CbmcGrow, F, GenCanPack, Guarantees, Handler, InsideBoxRestraint, Invariant,
    LatticeGrow, Layers, OnViolation, PackContext, PackEngine, PackError, PackSettings, PackState,
    Pipeline, Placed, Requires, RestraintsSatisfied, Stage, StageFactory, StageOutcome, State,
    StepInfo, Target, Until, Violation,
};
use molrs::store::block::Block;
use molrs::store::frame::Frame;
use ndarray::Array1;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

// ── templates (copied from tests/grow.rs and tests/packer.rs) ─────────────

/// Planar zigzag bead-chain coordinates: tetrahedral (109.5°) bond angles in
/// the x–z plane, all torsions trans. Copied from `tests/grow.rs` so this
/// file's goldens are independent of that file.
fn zigzag_coords(n: usize, bond_len: F) -> Vec<[F; 3]> {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let alpha = (std::f64::consts::PI as F - theta) / 2.0;
    let (dx, dz) = (bond_len * alpha.cos(), bond_len * alpha.sin());
    (0..n)
        .map(|i| [i as F * dx, 0.0, if i % 2 == 0 { 0.0 } else { dz }])
        .collect()
}

/// Zigzag bead chain with its `(i, i+1)` bond list, as a `molrs::Frame`
/// (atoms block with x/y/z, bonds block with atomi/atomj — the shape a PDB
/// CONECT list or a hand-built coarse-grain frame has).
fn chain_frame(n: usize, bond_len: F) -> Frame {
    let coords = zigzag_coords(n, bond_len);
    let mut atoms = Block::new();
    for (name, k) in [("x", 0), ("y", 1), ("z", 2)] {
        let col: Vec<F> = coords.iter().map(|p| p[k]).collect();
        atoms
            .insert(name, Array1::from_vec(col).into_dyn())
            .expect("coordinate column");
    }
    let mut frame = Frame::new();
    frame.insert("atoms", atoms);

    let mut bonds = Block::new();
    let ai: Vec<u32> = (0..n as u32 - 1).collect();
    let aj: Vec<u32> = (1..n as u32).collect();
    bonds
        .insert("atomi", Array1::from_vec(ai).into_dyn())
        .expect("atomi column");
    bonds
        .insert("atomj", Array1::from_vec(aj).into_dyn())
        .expect("atomj column");
    frame.insert("bonds", bonds);
    frame
}

/// Rigid water template (positions + packing radii), copied from
/// `tests/packer.rs`.
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
            .with_restraint(InsideBoxRestraint::new([0.0; 3], [14.0; 3], [false; 3])),
    ]
}

fn boxfree_settings<E: PackEngine>(engine: E) -> E {
    engine.with_seed(FREE_SEED).with_tolerance(FREE_TOL)
}

/// The growth fixture, mirroring `tests/grow.rs::seeded_run_contract` so the
/// `CbmcGrow` → `GenCanPack::with_restart` comparison is known to be
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

/// The melt fixture the abort tests use: eight 12-bead chains in a 26 Å
/// periodic box, template bond 1.53 Å (`tests/grow.rs`'s abort fixture).
const MELT_BOX: F = 26.0;
const MELT_BOND: F = 1.53;
const MELT_BEADS: usize = 12;
const MELT_COPIES: usize = 8;
const MELT_SEED: u64 = 7;
const MELT_TOL: F = 2.0;
const MELT_LOOPS: usize = 50;

fn melt_targets() -> Vec<Target> {
    vec![Target::new(chain_frame(MELT_BEADS, MELT_BOND), MELT_COPIES)]
}

fn melt_settings<E: PackEngine>(engine: E) -> E {
    engine
        .with_seed(MELT_SEED)
        .with_tolerance(MELT_TOL)
        .with_periodic_box([0.0; 3], [MELT_BOX; 3], [true; 3])
}

/// The lattice fixture: four 6-bead chains decorated onto a diamond lattice
/// in a 20 Å periodic box.
const LATTICE_SEED: u64 = 7;
const LATTICE_TOL: F = 2.0;
const LATTICE_LOOPS: usize = 40;

fn lattice_targets() -> Vec<Target> {
    vec![Target::new(chain_frame(6, 1.53), 4)]
}

fn lattice_settings<E: PackEngine>(engine: E) -> E {
    engine
        .with_seed(LATTICE_SEED)
        .with_tolerance(LATTICE_TOL)
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

/// Everything a run told its handlers, in arrival order.
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
    /// `StepInfo.stage` of every `on_step`.
    steps: Vec<(usize, usize, &'static str)>,
}

type Shared = Arc<Mutex<Tally>>;

fn shared() -> Shared {
    Arc::new(Mutex::new(Tally::default()))
}

/// Records every lifecycle callback into a shared [`Tally`], and — when
/// `stop_after` is finite — asks the run to stop once that many `on_step`
/// events have arrived (the `EarlyStopHandler` shape of
/// `tests/grow.rs::Recorder`).
struct Observer {
    tally: Shared,
    stop_after: usize,
}

impl Handler for Observer {
    fn on_start(&mut self, _ntotat: usize, _ntotmol: usize) {
        self.tally.lock().expect("observer mutex").starts += 1;
    }

    fn on_step(&mut self, info: &StepInfo, _sys: &PackContext) {
        self.tally.lock().expect("observer mutex").steps.push((
            info.stage.index,
            info.stage.total,
            info.stage.name,
        ));
    }

    fn on_stage_start(&mut self, info: &StageInfo) {
        self.tally
            .lock()
            .expect("observer mutex")
            .stage_starts
            .push((info.index, info.total, info.name));
    }

    fn on_stage_end(&mut self, info: &StageInfo, _outcome: &StageOutcome, _sys: &PackContext) {
        self.tally
            .lock()
            .expect("observer mutex")
            .stage_ends
            .push((info.index, info.total, info.name));
    }

    fn on_finish(&mut self, _sys: &PackContext) {
        self.tally.lock().expect("observer mutex").finishes += 1;
    }

    fn should_stop(&self) -> bool {
        self.tally.lock().expect("observer mutex").steps.len() >= self.stop_after
    }
}

/// An observer that never asks for a stop.
fn observer(tally: &Shared) -> Box<dyn Handler> {
    Box::new(Observer {
        tally: Arc::clone(tally),
        stop_after: usize::MAX,
    })
}

/// An observer that asks for a stop once `after` `on_step` events arrived.
fn stopping_observer(tally: &Shared, after: usize) -> Box<dyn Handler> {
    Box::new(Observer {
        tally: Arc::clone(tally),
        stop_after: after,
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
        _handlers: &mut [Box<dyn Handler>],
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

/// A single-stage GENCAN pipeline and the `GenCanPack` preset run are the
/// same arithmetic, bit for bit (ac-006). Declared periodic box.
#[test]
fn pipeline_single_stage_gencan_matches_preset_bitwise() {
    let targets = dense_targets();
    let piped = dense_settings(Pipeline::new().with_stage(GenCanPack::new()))
        .run(&targets, DENSE_LOOPS)
        .expect("a single-stage GENCAN pipeline runs");
    let direct = dense_settings(GenCanPack::new())
        .run(&targets, DENSE_LOOPS)
        .expect("the GenCanPack preset runs");

    assert_bitwise_equal(&piped, &direct, "single-stage GENCAN");
    assert_eq!(
        piped.converged, direct.converged,
        "the two spellings must also agree on the verdict"
    );
    assert_eq!(
        piped.softened, direct.softened,
        "a one-stage pipeline sums exactly one stage's softening"
    );
}

/// The same claim for the continuum growth entry (ac-006).
#[test]
fn pipeline_single_stage_cbmc_matches_preset_bitwise() {
    let targets = chain_targets();
    let piped = chain_settings(Pipeline::new().with_stage(CbmcGrow::new(TorsionPrior::Uniform)))
        .run(&targets, CHAIN_LOOPS)
        .expect("a single-stage CBMC pipeline runs");
    let direct = chain_settings(CbmcGrow::new(TorsionPrior::Uniform))
        .run(&targets, CHAIN_LOOPS)
        .expect("the CbmcGrow preset runs");

    assert_bitwise_equal(&piped, &direct, "single-stage CBMC growth");
    assert_eq!(piped.converged, direct.converged);
    assert_eq!(
        piped.softened, direct.softened,
        "the pipeline's softened total is the growth stage's own count"
    );
}

/// The same claim for the diamond-lattice growth entry (ac-006).
#[test]
fn pipeline_single_stage_lattice_matches_preset_bitwise() {
    let targets = lattice_targets();
    let piped =
        lattice_settings(Pipeline::new().with_stage(LatticeGrow::new(TorsionPrior::Uniform)))
            .run(&targets, LATTICE_LOOPS)
            .expect("a single-stage lattice pipeline runs");
    let direct = lattice_settings(LatticeGrow::new(TorsionPrior::Uniform))
        .run(&targets, LATTICE_LOOPS)
        .expect("the LatticeGrow preset runs");

    assert_bitwise_equal(&piped, &direct, "single-stage lattice growth");
    assert_eq!(piped.converged, direct.converged);
    assert_eq!(piped.softened, direct.softened);
}

// ── 2. cross-algorithm hand-off ≡ with_restart ─────────────────────────────

/// `[CbmcGrow, GenCanPack]` in a pipeline is bitwise
/// `GenCanPack::with_restart(&cbmc_result)` (ac-007).
///
/// The grown chains are the state the GENCAN stage inherits: it must
/// continue from them (push-off), never call `initial()` and scatter the
/// conformers growth just paid for. The seeded spelling is the same run
/// written by hand, so equality here is the whole point — and it only holds
/// if the pipeline invalidates the geometry cache at the stage boundary, so
/// both spellings enter GENCAN cold.
#[test]
fn pipeline_cbmc_then_gencan_equals_with_restart_bitwise() {
    let targets = chain_targets();

    let grown = chain_settings(CbmcGrow::new(TorsionPrior::Uniform))
        .run(&targets, CHAIN_LOOPS)
        .expect("the growth stage runs");
    assert!(
        grown.converged,
        "the seeded_run_contract fixture must grow cleanly for the hand-off \
         comparison to mean anything"
    );

    let piped = chain_settings(
        Pipeline::new()
            .with_stage(CbmcGrow::new(TorsionPrior::Uniform))
            .with_stage(GenCanPack::new()),
    )
    .run(&targets, CHAIN_LOOPS)
    .expect("a growth → GENCAN pipeline runs");

    // The cell travels with the seed, so the seeded spelling declares no box.
    let seeded = GenCanPack::new()
        .with_restart(&grown)
        .with_seed(CHAIN_SEED)
        .with_tolerance(CHAIN_TOL)
        .run(&targets, CHAIN_LOOPS)
        .expect("the hand-written seeded spelling runs");

    assert_bitwise_equal(&piped, &seeded, "[CbmcGrow, GenCanPack] vs with_restart");
}

/// `[GenCanPack, GenCanPack]` on one small budget is bitwise
/// `GenCanPack::with_restart(&first_result)` (ac-008).
///
/// The first stage is guarded to be UNCONVERGED, which is what makes the
/// claim observable: if the second GENCAN stage re-ran `initial()` it would
/// throw away the first stage's placements and land somewhere else entirely.
#[test]
fn pipeline_gencan_then_gencan_equals_with_restart_bitwise() {
    let targets = dense_targets();

    let first = dense_settings(GenCanPack::new())
        .run(&targets, DENSE_LOOPS)
        .expect("the first GENCAN stage runs");
    assert!(
        !first.converged,
        "guard: this fixture must NOT converge in {DENSE_LOOPS} loop(s) \
         (fdist = {}), otherwise a second stage that re-initialised would be \
         indistinguishable from one that continued",
        first.fdist
    );

    let piped = dense_settings(
        Pipeline::new()
            .with_stage(GenCanPack::new())
            .with_stage(GenCanPack::new()),
    )
    .run(&targets, DENSE_LOOPS)
    .expect("a GENCAN → GENCAN pipeline runs");

    let seeded = GenCanPack::new()
        .with_restart(&first)
        .with_seed(DENSE_SEED)
        .with_tolerance(DENSE_TOL)
        .run(&targets, DENSE_LOOPS)
        .expect("the hand-written seeded spelling runs");

    assert_bitwise_equal(
        &piped,
        &seeded,
        "[GenCanPack, GenCanPack] vs with_restart(first)",
    );
}

// ── 3. handlers are adopted, never dropped ────────────────────────────────

/// A handler attached to a PRESET that is then handed to `with_stage` still
/// receives its callbacks: the pipeline adopts it (ac-005). Dropping it is
/// the silent failure this whole design exists to prevent.
#[test]
fn pipeline_adopts_preset_handlers() {
    let one = shared();
    let piped = boxfree_settings(
        Pipeline::new().with_stage(GenCanPack::new().with_handler(observer(&one))),
    )
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a single-stage pipeline with an adopted handler runs");
    assert!(piped.converged);
    {
        let t = one.lock().expect("observer mutex");
        assert!(
            !t.steps.is_empty(),
            "the handler carried in by GenCanPack::with_handler saw no \
             on_step events — with_stage must adopt a preset's handlers, not \
             drop them on the floor"
        );
        assert_eq!(
            t.starts, 1,
            "an adopted handler is bracketed by the run like any other"
        );
        assert_eq!(t.finishes, 1);
    }

    // Two stages: an adopted handler observes the WHOLE run, so both stage
    // indices show up on the events it recorded.
    let both = shared();
    dense_settings(
        Pipeline::new()
            .with_stage(GenCanPack::new().with_handler(observer(&both)))
            .with_stage(GenCanPack::new()),
    )
    .run(&dense_targets(), DENSE_LOOPS)
    .expect("a two-stage pipeline with an adopted handler runs");

    let t = both.lock().expect("observer mutex");
    let mut seen: Vec<usize> = t.steps.iter().map(|&(index, _, _)| index).collect();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(
        seen,
        vec![0, 1],
        "an adopted handler must see BOTH stages (recorded stage indices \
         {seen:?}) — adoption is for the run, not for the stage that carried \
         it in"
    );
}

// ── 4. the two named rejections and the empty pipeline ────────────────────

/// An unmet `requires()` is reported by name BEFORE any handler is notified
/// and before any stage runs (ac-003).
#[test]
fn pipeline_stage_order_error_fires_before_any_handler() {
    let tally = shared();
    let err = boxfree_settings(
        Pipeline::new()
            .with_handler(observer(&tally))
            .with_stage(NeedsPlacedFactory::new())
            .with_stage(GenCanPack::new()),
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
         before ANY handler is notified"
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
        .with_stage(GenCanPack::new().with_seed(7))
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
    let direct = GenCanPack::new()
        .with_seed(7)
        .with_tolerance(FREE_TOL)
        .run(&boxfree_targets(), 1)
        .expect("a preset's own run adopts its settings and must not error");
    assert!(direct.fdist.is_finite() && direct.frest.is_finite());

    // …and neither is the explicit `Pipeline::single` spelling.
    let single = Pipeline::single(GenCanPack::new().with_seed(7).with_tolerance(FREE_TOL))
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

// ── 5. one verdict, one bracket per run, one bracket per stage ────────────

/// A two-stage run brackets the RUN once and each STAGE once, reports
/// `stage.total == 2` on every step with a non-decreasing `stage.index`, and
/// sums `softened` across the stages (ac-009).
#[test]
fn pipeline_two_stages_sum_softened_and_count_hooks() {
    let targets = chain_targets();

    // The two standalone spellings, same seed and settings, provide the
    // reference sum.
    let grown = chain_settings(CbmcGrow::new(TorsionPrior::Uniform))
        .run(&targets, CHAIN_LOOPS)
        .expect("the growth stage runs standalone");
    let seeded = GenCanPack::new()
        .with_restart(&grown)
        .with_seed(CHAIN_SEED)
        .with_tolerance(CHAIN_TOL)
        .run(&targets, CHAIN_LOOPS)
        .expect("the seeded GENCAN stage runs standalone");

    let carried = shared();
    let attached = shared();
    let piped = chain_settings(
        Pipeline::new()
            .with_handler(observer(&attached))
            .with_stage(CbmcGrow::new(TorsionPrior::Uniform).with_handler(observer(&carried)))
            .with_stage(GenCanPack::new()),
    )
    .run(&targets, CHAIN_LOOPS)
    .expect("a growth → GENCAN pipeline runs");

    assert_eq!(
        piped.softened,
        grown.softened + seeded.softened,
        "the pipeline's softened count is the SUM over its stages \
         ({} + {}), not the last stage's own count",
        grown.softened,
        seeded.softened
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
            "{label}: each StageInfo carries the stage's own name()"
        );
        for &(index, total, _) in t.stage_starts.iter().chain(t.stage_ends.iter()) {
            assert_eq!(
                total, 2,
                "{label}: StageInfo.total is the number of stages in the run"
            );
            assert!(index < total, "{label}: stage index {index} out of range");
        }
        for &(index, total, _) in t.steps.iter() {
            assert_eq!(
                total, 2,
                "{label}: every StepInfo reports the pipeline's stage count"
            );
            assert!(index < total);
        }
        for pair in t.steps.windows(2) {
            let (prev, cur) = (pair[0].0, pair[1].0);
            assert!(
                cur >= prev,
                "{label}: StepInfo.stage.index went {prev} → {cur} — the \
                 stage index is monotone along a linear pipeline"
            );
        }
    }
}

/// A `should_stop` inside the FIRST of two stages ends the whole run: the
/// second stage never starts, the result is honestly unconverged, and the
/// growth stage's abort contract still leaves chemical 1-2 distances
/// (ac-009, the predicate of `grow_abort_keeps_bonded_geometry`).
#[test]
fn pipeline_early_stop_keeps_bonded_geometry_and_is_unconverged() {
    let tally = shared();
    let piped = melt_settings(
        Pipeline::new()
            .with_handler(stopping_observer(&tally, 2))
            .with_stage(CbmcGrow::new(TorsionPrior::Uniform))
            .with_stage(GenCanPack::new()),
    )
    .run(&melt_targets(), MELT_LOOPS)
    .expect("an aborted pipeline still returns Ok");

    assert!(
        !piped.converged,
        "a run abandoned mid-way must never be reported as converged"
    );
    {
        let t = tally.lock().expect("observer mutex");
        assert_eq!(
            t.stage_starts.len(),
            1,
            "should_stop in stage 0 must end the run: the second stage's \
             on_stage_start fired anyway ({:?})",
            t.stage_starts
        );
        assert_eq!(
            t.stage_ends.len(),
            1,
            "the aborted stage is still closed exactly once"
        );
        assert_eq!(t.starts, 1, "the run bracket still opened once");
        assert_eq!(
            t.finishes, 1,
            "an aborted run still closes its handler bracket"
        );
    }

    // Bonded geometry survives the abort: an unplaced atom left at the
    // origin would show up as a 0-length or box-scale 1-2 distance.
    let pos = piped.positions();
    assert_eq!(pos.len(), MELT_COPIES * MELT_BEADS);
    let bonds = piped
        .frame
        .get("bonds")
        .expect("the assembled frame keeps the tiled template bonds");
    let ai = bonds.get_uint("atomi").expect("atomi column");
    let aj = bonds.get_uint("atomj").expect("atomj column");
    assert!(
        !ai.is_empty(),
        "template 1-2 bonds must be tiled onto copies"
    );
    for (&a, &b) in ai.iter().zip(aj.iter()) {
        let (p, q) = (pos[a as usize], pos[b as usize]);
        let d: [F; 3] = std::array::from_fn(|k| {
            let raw = p[k] - q[k];
            raw - (raw / MELT_BOX).round() * MELT_BOX
        });
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        assert!(
            (len - MELT_BOND).abs() < 1e-9,
            "1-2 distance {a}–{b} is {len} Å, not the template bond \
             {MELT_BOND} Å — an aborted growth stage must complete its \
             chains chemically, never leave origin sentinels behind"
        );
    }
}

// ── 6. the box-free single-stage run ──────────────────────────────────────

/// A single-stage GENCAN pipeline with NO box and NO cell declaration runs
/// to completion (ac-006): the GENCAN stage's preamble installs nothing
/// while the state is `Placed::None` and unseeded, leaving `initial()` to
/// synthesize its fall-back box — and it is still bitwise the preset run.
#[test]
fn pipeline_gencan_without_box_runs() {
    let targets = boxfree_targets();
    let piped = boxfree_settings(Pipeline::new().with_stage(GenCanPack::new()))
        .run(&targets, FREE_LOOPS)
        .expect("a box-free single-stage GENCAN pipeline must run, not panic");
    let direct = boxfree_settings(GenCanPack::new())
        .run(&targets, FREE_LOOPS)
        .expect("the box-free preset run");

    assert!(
        piped.fdist.is_finite() && piped.frest.is_finite(),
        "box-free run produced non-finite violations: fdist = {}, frest = {}",
        piped.fdist,
        piped.frest
    );
    for (i, p) in piped.positions().iter().enumerate() {
        for (k, x) in p.iter().enumerate() {
            assert!(x.is_finite(), "atom {i} component {k} is {x}");
        }
    }
    assert_bitwise_equal(&piped, &direct, "box-free single-stage GENCAN");
}

// ── 7. regression scenario — hard-coded goldens ───────────────────────────

/// Hard-coded golden for the box-free single-stage GENCAN pipeline
/// (ac-012).
///
/// Provenance: captured 2026-09-03 from the build at commit 77cba83, before
/// stage-pipeline-05. Tool: this repository's own `cargo test
/// -p molcrafts-molpack` (debug profile, default features) running a scratch
/// integration test that printed `State::fdist` / `frest` and
/// `positions()` through `{:?}` (shortest round-trip form) for exactly the
/// fixture below, spelled as `GenCanPack::new().with_seed(42)
/// .with_tolerance(2.0).run(&boxfree_targets(), 20)`. That preset spelling
/// is the pipeline spelling by ac-006, which is why a pre-pipeline capture
/// is the right oracle. No third-party program was involved.
///
/// Deterministic by construction: fixed seed 42, serial evaluation, no wall
/// clock, no filesystem, no network.
#[test]
fn pipeline_regression_single_stage_gencan_golden() {
    /// Largest inter-molecular contact violation at termination — a strict
    /// zero on this generously feasible fixture.
    const GOLDEN_FDIST: F = 0.0;
    /// Largest restraint violation at termination.
    const GOLDEN_FREST: F = 0.000_301_283_868_137_773_27;
    /// The first water's O, H, H in target-declared order (Å).
    const GOLDEN_HEAD: [[F; 3]; 3] = [
        [13.068219535648394, 9.39835502863678, 10.821668445710404],
        [13.389024780022366, 8.517418653690637, 11.028151128952619],
        [12.81225014697141, 9.765834204264351, 11.671338221294079],
    ];
    const TOL: F = 1e-12;

    let result = boxfree_settings(Pipeline::new().with_stage(GenCanPack::new()))
        .run(&boxfree_targets(), FREE_LOOPS)
        .expect("the golden fixture runs");

    assert!(
        (result.fdist - GOLDEN_FDIST).abs() < TOL,
        "fdist drifted: {} vs golden {GOLDEN_FDIST}",
        result.fdist
    );
    assert!(
        (result.frest - GOLDEN_FREST).abs() < TOL,
        "frest drifted: {} vs golden {GOLDEN_FREST}",
        result.frest
    );

    let pos = result.positions();
    assert_eq!(pos.len(), 180, "60 waters × 3 atoms");
    for (i, (got, want)) in pos.iter().take(3).zip(GOLDEN_HEAD.iter()).enumerate() {
        for (k, (g, w)) in got.iter().zip(want.iter()).enumerate() {
            assert!(
                (g - w).abs() < TOL,
                "atom {i} component {k}: {g} vs golden {w} (|Δ| ≥ {TOL})"
            );
        }
    }
}

// ── 8. the two combinators (spec 06) ──────────────────────────────────────
//
// `Repeat` and `Guarded` are stages, so everything above still applies to
// them: they are resolved by the same chain check, bracketed by the same
// handlers, and report ONE stage identity no matter how often their body
// runs. What is new here is the body's run count, the honest `softened` sum,
// the continuation of a repeated GENCAN pass, and the named refusal a broken
// invariant produces.

/// A stage that moves no atoms and counts its runs through a shared counter.
///
/// The combinators are about *how often* and *under what condition* a stage
/// runs, so the body they wrap only has to be countable — a real algorithm
/// would add arithmetic that hides the count. When asked, the stage also
/// emits one `on_inner_iter` per run: that is how a handler learns the body
/// did a unit of work from *inside* a combinator, because
/// [`StepInfo`](molpack::StepInfo) is `#[non_exhaustive]` and an integration
/// test therefore cannot build one to drive `on_step`.
struct CountingStage {
    runs: Arc<AtomicUsize>,
    converged: bool,
    softened: usize,
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
        state: &mut PackState,
        _targets: &[Target],
        _budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        self.runs.fetch_add(1, Ordering::Relaxed);
        if self.signal {
            for h in handlers.iter_mut() {
                h.on_inner_iter(0, 0.0, state.ctx());
            }
        }
        Ok(StageOutcome::new(self.converged, self.softened))
    }
}

/// The factory that produces one [`CountingStage`], carrying default shared
/// settings so a pipeline has nothing to refuse.
struct CountingFactory {
    settings: PackSettings,
    runs: Arc<AtomicUsize>,
    converged: bool,
    softened: usize,
    signal: bool,
}

impl CountingFactory {
    fn new(runs: &Arc<AtomicUsize>) -> Self {
        Self {
            settings: PackSettings::default(),
            runs: Arc::clone(runs),
            converged: false,
            softened: 0,
            signal: false,
        }
    }

    /// The produced stage reports its own convergence criterion as met.
    fn converging(mut self) -> Self {
        self.converged = true;
        self
    }

    /// The produced stage reports `n` relaxations per run.
    fn with_softened(mut self, n: usize) -> Self {
        self.softened = n;
        self
    }

    /// The produced stage emits one `on_inner_iter` per run.
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
            softened: self.softened,
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

/// A handler that asks for a stop as soon as the body has signalled one unit
/// of work. `should_stop` must be honoured from *inside* a combinator, not
/// only between the pipeline's own stages.
struct StopOnFirstSignal {
    signals: Arc<AtomicUsize>,
}

impl Handler for StopOnFirstSignal {
    fn on_step(&mut self, _info: &StepInfo, _sys: &PackContext) {}

    fn on_inner_iter(&mut self, _iter: u32, _f: F, _sys: &PackContext) {
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
fn repeat_passes_runs_body_n_times_and_sums_softened() {
    const SOFTENED_PER_RUN: usize = 3;
    let runs = Arc::new(AtomicUsize::new(0));

    let result = boxfree_settings(Pipeline::new().with_repeat(
        body(CountingFactory::new(&runs).with_softened(SOFTENED_PER_RUN)),
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
        result.softened,
        2 * SOFTENED_PER_RUN,
        "the repeated stage's softened count is the SUM over its passes \
         (2 × {SOFTENED_PER_RUN}), not one pass's own count"
    );
}

/// The core of ac-002: the second pass **continues** from the first.
///
/// `Repeat{Passes(2)}` around one GENCAN stage is bitwise
/// `GenCanPack::with_restart(&pass1)` — the same comparison
/// `pipeline_gencan_then_gencan_equals_with_restart_bitwise` makes for two
/// explicit stages. Without it, "repeat" and "run the whole thing again from
/// a fresh `initial()`" are indistinguishable in a test.
///
/// **Why bitwise equality is the whole proof, and `fdist` is not part of
/// it.** One more outer loop is *not* monotone in the unscaled `fdist`:
/// `run_phase` restarts `radscale` at `discale` for every phase, so a
/// continuation legitimately climbs before it descends again. Measured on
/// this fixture, pass 1 ends at `fdist = 9.0308` and all three continuation
/// spellings — `with_restart`, `[GenCanPack, GenCanPack]` and this `Repeat` —
/// end at `9.7624`, bit for bit. Agreeing with the seeded spelling to the
/// bit is therefore the assertion that separates "continued" from
/// "re-`initial()`ed"; a `fdist` inequality would only be asserting the
/// radius schedule.
#[test]
fn repeat_second_pass_continues_from_the_first_bitwise() {
    let targets = dense_targets();

    let pass1 = dense_settings(GenCanPack::new())
        .run(&targets, DENSE_LOOPS)
        .expect("the first GENCAN pass runs");
    assert!(
        !pass1.converged,
        "guard: this fixture must NOT converge in {DENSE_LOOPS} loop(s) \
         (fdist = {}), otherwise a second pass that re-initialised would be \
         indistinguishable from one that continued",
        pass1.fdist
    );

    let repeated =
        dense_settings(Pipeline::new().with_repeat(body(GenCanPack::new()), Until::Passes(2)))
            .run(&targets, DENSE_LOOPS)
            .expect("a Repeat{Passes(2)} GENCAN pipeline runs");

    let seeded = GenCanPack::new()
        .with_restart(&pass1)
        .with_seed(DENSE_SEED)
        .with_tolerance(DENSE_TOL)
        .run(&targets, DENSE_LOOPS)
        .expect("the hand-written seeded spelling runs");

    assert_bitwise_equal(
        &repeated,
        &seeded,
        "Repeat{Passes(2)} vs with_restart(pass1)",
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
            .with_stage(GenCanPack::new())
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

/// A guard that holds changes nothing: `Guarded(stage, [invariant], _)` with
/// a satisfied invariant is bitwise the bare stage. A combinator that
/// perturbed the run it observes would make every guarded result a different
/// run from the unguarded one.
#[test]
fn guarded_passing_invariant_leaves_result_bitwise() {
    let targets = boxfree_targets();
    // The box-free fixture ends at frest ≈ 3.0e-4 (pinned by
    // `pipeline_regression_single_stage_gencan_golden`), well inside 1.0.
    let invariants: Vec<Box<dyn Invariant>> = vec![Box::new(RestraintsSatisfied::new(1.0))];

    let guarded = boxfree_settings(Pipeline::new().with_guarded(
        GenCanPack::new(),
        invariants,
        OnViolation::Fail,
    ))
    .run(&targets, FREE_LOOPS)
    .expect("a guarded GENCAN stage whose invariant holds runs");
    let plain = boxfree_settings(Pipeline::new().with_stage(GenCanPack::new()))
        .run(&targets, FREE_LOOPS)
        .expect("the unguarded spelling runs");

    assert_bitwise_equal(&guarded, &plain, "Guarded(gencan, [satisfied]) vs gencan");
    assert_eq!(
        guarded.converged, plain.converged,
        "a satisfied guard must not change the verdict either"
    );
}

/// `OnViolation::Fail` surfaces a NAMED error out of `Pipeline::run`, naming
/// the stage, the invariant, its rung on the repair-cost ladder and the atoms
/// involved — never a quietly unconverged result and never another algorithm
/// (law P8).
///
/// The atom list itself is exercised on a hand-built state in
/// `tests/invariant.rs`: `frest_atom` is movebad-only bookkeeping
/// (`src/movebad.rs`), so which atoms a real run's final state attributes the
/// residual to is the GENCAN schedule's business, not this file's.
#[test]
fn guarded_fail_returns_named_error() {
    // frest ≈ 3.0e-4 on this fixture, so a tolerance of exactly 0.0 is
    // unsatisfiable by construction.
    let invariants: Vec<Box<dyn Invariant>> = vec![Box::new(RestraintsSatisfied::new(0.0))];

    let err = boxfree_settings(Pipeline::new().with_guarded(
        GenCanPack::new(),
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
/// move `StepInfo.stage.index` or `total`, and the names a handler sees are
/// the top-level ones. Nesting a combinator inside another does not change
/// that (ac-006).
#[test]
fn combinators_keep_stage_index_monotone() {
    let flat = shared();
    dense_settings(
        Pipeline::new()
            .with_handler(observer(&flat))
            .with_stage(GenCanPack::new())
            .with_repeat(body(GenCanPack::new()), Until::Passes(2)),
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
                "StepInfo.stage.index went {prev} → {cur} — a repeated body \
                 must not walk the index backwards"
            );
        }
    }

    // Nested: Guarded(Repeat([gencan])). The invariant list is empty on
    // purpose — this test is about stage identity, and the guard's own
    // verdict is pinned by the two tests above.
    let nested = shared();
    let inner = Pipeline::new().with_repeat(body(GenCanPack::new()), Until::Passes(2));
    dense_settings(
        Pipeline::new()
            .with_handler(observer(&nested))
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

/// A body factory's handlers are **adopted** exactly as `with_stage` adopts
/// them, and its non-default shared settings are refused BY NAME exactly as
/// `with_stage` refuses them — a combinator is not a hole in either rule.
#[test]
fn repeat_adopts_body_handlers_and_refuses_settings() {
    let tally = shared();
    boxfree_settings(Pipeline::new().with_repeat(
        body(GenCanPack::new().with_handler(observer(&tally))),
        Until::Passes(2),
    ))
    .run(&boxfree_targets(), FREE_LOOPS)
    .expect("a Repeat whose body carries a handler runs");

    {
        let t = tally.lock().expect("observer mutex");
        assert!(
            !t.steps.is_empty(),
            "the handler carried in by the body preset saw no on_step events \
             — with_repeat must adopt a body's handlers, not drop them"
        );
        assert_eq!(t.starts, 1, "an adopted handler is bracketed once per RUN");
        assert_eq!(t.finishes, 1);
    }

    let err = Pipeline::new()
        .with_repeat(body(GenCanPack::new().with_seed(7)), Until::Passes(2))
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

/// A handler asking to stop from inside a combinator body ends the run there:
/// the remaining passes never start and the verdict is honestly
/// `converged == false`.
#[test]
fn repeat_should_stop_yields_immediately() {
    let runs = Arc::new(AtomicUsize::new(0));
    let signals = Arc::new(AtomicUsize::new(0));

    let result = boxfree_settings(
        Pipeline::new()
            .with_handler(Box::new(StopOnFirstSignal {
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
