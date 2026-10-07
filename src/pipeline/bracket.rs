//! The callback bracket around a run.
//!
//! Adopted-callback tagging ([`StageTagger`], [`tag_all`]) stamps the run's
//! current stage position onto every `on_step` event a callback receives —
//! a stage only knows `index: 0, total: 1`, so where it sits in a chain has
//! to be stamped on from outside. [`open_bracket`] opens the bracket:
//! appends the built-in LAMMPS log callback when enabled, tags the set, and
//! fires `on_start`. [`close_bracket`] closes it with `on_finish`. The
//! per-stage loop that runs between the two is the lifecycle body and lives
//! in `mod.rs`, not here.

use std::sync::{Arc, Mutex};

use crate::callback::{
    Callback, LammpsLogCallback, PhaseProgress, PhaseReport, StageProgress, StepReport,
};
use crate::context::PackContext;
use crate::pack_space::ResolvedSpace;
use crate::stage::{Budget, Stage, StageOutcome};

use super::EngineSetup;

/// The slot the pipeline writes its current stage identity into, shared with
/// every wrapped callback.
pub(super) type StagePosition = Arc<Mutex<StageProgress>>;

/// Wrap each callback so its `on_step` events carry the run's stage position.
///
/// One pass over the set, at the point the run's callback list is assembled;
/// nothing else in the crate wraps a callback.
pub(super) fn tag_all(
    callbacks: Vec<Box<dyn Callback>>,
    position: &StagePosition,
) -> Vec<Box<dyn Callback>> {
    callbacks
        .into_iter()
        .map(|inner| {
            Box::new(StageTagger {
                inner,
                position: Arc::clone(position),
            }) as Box<dyn Callback>
        })
        .collect()
}

/// Open the run's callback bracket: append the built-in LAMMPS log callback
/// when enabled, tag the set with a fresh position slot (stamped at the
/// first stage), fire `on_start`, and hand both back —
/// [`Pipeline::run_stages`](super::Pipeline::run_stages) updates the slot
/// before each stage.
pub(super) fn open_bracket(
    own_callbacks: Vec<Box<dyn Callback>>,
    stages: &[Box<dyn Stage>],
    setup: &EngineSetup<'_>,
    space: &ResolvedSpace,
    budget: &Budget,
) -> (Vec<Box<dyn Callback>>, StagePosition) {
    let log = setup.settings.log;
    let mut collected = own_callbacks;
    if log.level.is_enabled() {
        collected.push(Box::new(LammpsLogCallback::new(
            log.level,
            log.frequency,
            setup.settings.tolerance(),
            budget.precision,
            setup.settings.seed(),
            budget.max_loops,
            setup.ntype_with_fixed,
            space.cell.clone(),
        )));
    }
    let position = Arc::new(Mutex::new(StageProgress {
        index: 0,
        total: stages.len(),
        name: stages[0].name(),
    }));
    let mut tagged = tag_all(collected, &position);
    for h in tagged.iter_mut() {
        h.on_start(setup.ntotat, setup.ntotmol_free);
    }
    (tagged, position)
}

/// Close the run's callback bracket: fire `on_finish` on the tagged set
/// [`open_bracket`] produced.
pub(super) fn close_bracket(callbacks: &mut [Box<dyn Callback>], sys: &PackContext) {
    for h in callbacks.iter_mut() {
        h.on_finish(sys);
    }
}

/// Stamps the pipeline's current stage position onto every `StepReport` a
/// stage emits.
///
/// A stage reports `index: 0, total: 1` — true standalone, and all it can
/// know. Where it sits in a chain is the pipeline's fact, so the pipeline
/// owns it: one wrapper per callback, one shared slot updated before each
/// stage runs. Every other callback, including `should_stop`, is forwarded
/// untouched.
struct StageTagger {
    inner: Box<dyn Callback>,
    position: StagePosition,
}

impl Callback for StageTagger {
    fn on_start(&mut self, ntotat: usize, ntotmol: usize) {
        self.inner.on_start(ntotat, ntotmol);
    }

    fn on_initialized(&mut self, sys: &PackContext) {
        self.inner.on_initialized(sys);
    }

    fn on_step(&mut self, step: &StepReport, sys: &PackContext) {
        let mut tagged: StepReport = step.clone();
        tagged.stage = *self.position.lock().expect("stage position mutex");
        self.inner.on_step(&tagged, sys);
    }

    fn on_phase_start(&mut self, phase: &PhaseProgress) {
        self.inner.on_phase_start(phase);
    }

    fn on_finish(&mut self, sys: &PackContext) {
        self.inner.on_finish(sys);
    }

    fn should_stop(&self) -> bool {
        self.inner.should_stop()
    }

    fn on_phase_end(&mut self, phase: &PhaseProgress, report: &PhaseReport) {
        self.inner.on_phase_end(phase, report);
    }

    fn on_stage_start(&mut self, stage: &StageProgress) {
        self.inner.on_stage_start(stage);
    }

    fn on_stage_end(&mut self, stage: &StageProgress, outcome: &StageOutcome, sys: &PackContext) {
        self.inner.on_stage_end(stage, outcome, sys);
    }
}
