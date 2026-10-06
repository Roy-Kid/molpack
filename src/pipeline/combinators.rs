//! Two ways to compose stages: repeat a body, or guard a stage's exit.
//!
//! **A combinator is a [`Stage`].** [`Pipeline`](super::Pipeline) needs no
//! branch for either — chain check, handler bracket and verdict read them as
//! they read `GenCanPack` — and each stays *one* stage to the run around it:
//! `name()` is `"repeat"` or `"guarded"`, inner passes never touch the index
//! and total the pipeline stamps on `StepInfo.stage`, and the inner stages
//! get no `on_stage_start` / `on_stage_end` of their own.
//!
//! What a combinator owes its body is the rest of `mod.rs`'s per-stage
//! boundary, copied from there rather than invented here (see `run_body`):
//! a cold geometry cache, and the placement marker advanced by each inner
//! stage's [`Guarantees`] — which is what makes pass *n+1* **continue** pass
//! *n* instead of packing again from nothing. Two earlier promises hold it
//! up: the re-entrancy contract on [`Stage::run`], and `Placed::All`
//! continuation. `requires()` is the body's first stage's, `Repeat`'s
//! `guarantees()` the union of the body's, `Guarded`'s the last stage's;
//! `degraded` sums over every pass and every attempt. Why the crate pays
//! for any of it, and how to roll it back: the ledger in [`mod`](super)'s
//! docs. What a guard may never do: [`OnViolation`].

use std::sync::OnceLock;

use crate::Handler;
use crate::Invariant;
use crate::PackError;
use crate::PackSettings;
use crate::Target;
use crate::context::{PackState, Placed};
use crate::stage::{Budget, Guarantees, Requires, Stage, StageOutcome};

use super::engine::{EngineSetup, StageFactory};

/// When [`Pipeline::with_repeat`](super::Pipeline::with_repeat) stops.
#[derive(Debug, Clone, Copy)]
pub enum Until {
    /// Exactly `n` passes. `Passes(0)` contributes **no stage at all** — the
    /// body never runs, and a pipeline left empty by it is the named
    /// [`PackError::NoStages`]; it is not clamped to one pass.
    Passes(usize),
    /// Until a pass ends with the body's last stage reporting its own
    /// convergence criterion met — unbounded by construction: a body that
    /// never converges repeats until a handler stops the run.
    Converged,
}

/// What [`Pipeline::with_guarded`](super::Pipeline::with_guarded) does when
/// an invariant is broken at the guarded stage's exit. Two arms, never a
/// third: a guard reruns **the same stage** or fails by name, because
/// picking the method is the caller's act (law P8).
#[derive(Debug, Clone, Copy)]
pub enum OnViolation {
    /// Fail the run with [`PackError::InvariantViolated`].
    Fail,
    /// Rerun the same stage at most `max` more times; still broken, report
    /// `converged == false` and log it. `Rerun { max: 0 }` is [`Fail`](Self::Fail).
    Rerun {
        /// How many *re*runs are allowed after the first attempt.
        max: usize,
    },
}

/// A combinator carries no shared knobs: the run's ruler is the pipeline's.
fn no_knobs() -> &'static PackSettings {
    static DEFAULTS: OnceLock<PackSettings> = OnceLock::new();
    DEFAULTS.get_or_init(PackSettings::default)
}

/// Resolve one body factory's stages, refusing a second ruler by name first
/// — a combinator is no hole in the one-ruler rule. `fallback` names the
/// combinator when the offending factory produced no stage to name.
fn resolve(
    factory: &mut Box<dyn StageFactory>,
    setup: &EngineSetup<'_>,
    fallback: &'static str,
) -> Result<Vec<Box<dyn Stage>>, PackError> {
    let produced = factory.stages(setup)?;
    match factory.settings().first_non_default_knob() {
        Some(knob) => Err(PackError::PresetSettingsInsidePipeline {
            stage: produced.first().map_or(fallback, |s| s.name()),
            knob,
        }),
        None => Ok(produced),
    }
}

/// Run `stages` once, in order, across the same boundary the pipeline puts
/// in front of a stage. Reports whether a handler asked to stop, which a
/// combinator honours immediately.
fn run_body(
    stages: &mut [Box<dyn Stage>],
    state: &mut PackState,
    targets: &[Target],
    budget: &Budget,
    handlers: &mut [Box<dyn Handler>],
) -> Result<(StageOutcome, bool), PackError> {
    let (mut converged, mut degraded) = (false, 0usize);
    for stage in stages.iter_mut() {
        state.invalidate_geometry_cache();
        let outcome = stage.run(state, targets, budget, handlers)?;
        state.set_placed(stage.guarantees().placed);
        degraded += outcome.degraded;
        converged = outcome.converged;
        if handlers.iter().any(|h| h.should_stop()) {
            return Ok((StageOutcome::new(false, degraded), true));
        }
    }
    Ok((StageOutcome::new(converged, degraded), false))
}

const REPEAT: &str = "repeat";
const GUARDED: &str = "guarded";

/// What [`Pipeline::with_repeat`](super::Pipeline::with_repeat) stores.
pub(super) struct RepeatFactory {
    pub(super) body: Vec<Box<dyn StageFactory>>,
    pub(super) until: Until,
}

impl StageFactory for RepeatFactory {
    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        self.body
            .iter()
            .try_for_each(|f| f.validate_targets(targets))
    }

    fn settings(&self) -> &PackSettings {
        no_knobs()
    }

    /// Adopted, in body order, exactly as `with_stage` adopts a preset's.
    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        let taken = self.body.iter_mut().map(|f| f.take_handlers());
        taken.flatten().collect()
    }

    /// One [`Repeat`], or none. The body is resolved (and refused) even for
    /// `Passes(0)`: one input, one verdict, whatever the stopping rule says.
    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        let mut stages: Vec<Box<dyn Stage>> = Vec::new();
        for factory in self.body.iter_mut() {
            stages.extend(resolve(factory, setup, REPEAT)?);
        }
        if stages.is_empty() || matches!(self.until, Until::Passes(0)) {
            return Ok(Vec::new());
        }
        Ok(vec![Box::new(Repeat {
            stages,
            until: self.until,
        })])
    }
}

/// The body, run until [`Until`] is met. `stages` is non-empty.
struct Repeat {
    stages: Vec<Box<dyn Stage>>,
    until: Until,
}

impl Stage for Repeat {
    fn name(&self) -> &'static str {
        REPEAT
    }

    fn requires(&self) -> Requires {
        self.stages[0].requires()
    }

    fn guarantees(&self) -> Guarantees {
        let mut body = self.stages.iter();
        let all = body.any(|s| s.guarantees().placed == Placed::All);
        Guarantees::new(if all { Placed::All } else { Placed::None })
    }

    fn run(
        &mut self,
        state: &mut PackState,
        targets: &[Target],
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        let (mut degraded, mut passes) = (0usize, 0usize);
        let mut converged;
        loop {
            let (outcome, stopped) = run_body(&mut self.stages, state, targets, budget, handlers)?;
            degraded += outcome.degraded;
            converged = outcome.converged;
            passes += 1;
            if stopped {
                return Ok(StageOutcome::new(false, degraded));
            }
            match self.until {
                Until::Converged if converged => break,
                Until::Passes(n) if passes >= n => break,
                _ => {}
            }
        }
        Ok(StageOutcome::new(converged, degraded))
    }
}

/// What [`Pipeline::with_guarded`](super::Pipeline::with_guarded) stores.
pub(super) struct GuardedFactory {
    pub(super) inner: Box<dyn StageFactory>,
    pub(super) invariants: Vec<Box<dyn Invariant>>,
    pub(super) on_violation: OnViolation,
}

impl StageFactory for GuardedFactory {
    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        self.inner.validate_targets(targets)
    }

    fn settings(&self) -> &PackSettings {
        no_knobs()
    }

    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        self.inner.take_handlers()
    }

    /// One [`Guarded`] around whatever the inner factory produced — a
    /// nested `Pipeline` may hand over several, guarded as one sequence.
    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        let stages = resolve(&mut self.inner, setup, GUARDED)?;
        if stages.is_empty() {
            return Ok(Vec::new());
        }
        Ok(vec![Box::new(Guarded {
            stages,
            invariants: std::mem::take(&mut self.invariants),
            on_violation: self.on_violation,
        })])
    }
}

/// The guarded stage(s) plus the invariants their exit must satisfy.
struct Guarded {
    stages: Vec<Box<dyn Stage>>,
    invariants: Vec<Box<dyn Invariant>>,
    on_violation: OnViolation,
}

impl Stage for Guarded {
    fn name(&self) -> &'static str {
        GUARDED
    }

    fn requires(&self) -> Requires {
        self.stages[0].requires()
    }

    fn guarantees(&self) -> Guarantees {
        self.stages[self.stages.len() - 1].guarantees()
    }

    /// Run the guarded stage(s), check the invariants on the state they
    /// left, answer by [`OnViolation`]. First broken invariant only.
    fn run(
        &mut self,
        state: &mut PackState,
        targets: &[Target],
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        let (mut degraded, mut attempt) = (0usize, 0usize);
        let stage = self.stages[0].name();
        loop {
            let (outcome, stopped) = run_body(&mut self.stages, state, targets, budget, handlers)?;
            degraded += outcome.degraded;
            if stopped {
                return Ok(StageOutcome::new(false, degraded));
            }
            let broken = self.invariants.iter().find_map(|inv| {
                let violation = inv.check(state).into_iter().next()?;
                Some((inv.name(), inv.layer().name(), violation.atoms))
            });
            let Some((invariant, layer, atoms)) = broken else {
                return Ok(StageOutcome::new(outcome.converged, degraded));
            };
            match self.on_violation {
                OnViolation::Rerun { max } if attempt < max => {
                    attempt += 1;
                    log::warn!("  `{stage}` broke `{invariant}` ({layer}); rerun {attempt}/{max}");
                }
                OnViolation::Fail | OnViolation::Rerun { max: 0 } => {
                    return Err(PackError::InvariantViolated {
                        stage,
                        invariant,
                        layer,
                        atoms,
                    });
                }
                OnViolation::Rerun { max } => {
                    log::warn!(
                        "  `{stage}` still breaks `{invariant}` ({layer}) after {max} rerun(s) \
                         — reporting converged = false; molpack does not switch algorithm"
                    );
                    return Ok(StageOutcome::new(false, degraded));
                }
            }
        }
    }
}
