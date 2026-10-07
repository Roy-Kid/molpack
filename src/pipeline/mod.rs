//! The run lifecycle: [`Pipeline`](crate::Pipeline) and its documentation.

mod bracket;
mod combinators;
mod engine;
#[cfg(test)]
mod tests;

use molrs::op::F;

use crate::Invariant;
use crate::PackError;
use crate::Target;
use crate::context::build::{ContextKnobs, build_context};
use crate::context::{PackState, Placed};
use crate::entry::setup::{ResolvedSpace, broadcast_global_restraints, resolve_pack_space};
use crate::entry::{IntraResidual, PackSettings, State, positions_in_target_order, result};
use crate::handler::{Handler, StageInfo};
use crate::stage::{Budget, Stage};
use bracket::{close_bracket, open_bracket};
use combinators::{GuardedFactory, RepeatFactory};

pub use combinators::{OnViolation, Until};
pub use engine::{EngineSetup, PackEngine, StageFactory};

/// A sequence of stages behind one lifecycle, one settings set and one
/// handler set.
///
/// ```no_run
/// use molpack::grow::TorsionPrior;
/// use molpack::{CbmcGrow, GenCanPack, PackEngine, Pipeline, Target};
/// # let targets: Vec<Target> = Vec::new();
/// let result = Pipeline::new()
///     .with_stage(CbmcGrow::new(TorsionPrior::Uniform))
///     .with_stage(GenCanPack::new())
///     .with_seed(42)
///     .with_periodic_box([0.0; 3], [30.0; 3], [true; 3])
///     .run(&targets, 100)?;
/// # Ok::<(), molpack::PackError>(())
/// ```
///
/// The lifecycle body: one run, one place.
///
/// [`Pipeline`] holds the shared settings, the run's handlers and a sequence
/// of [`StageFactory`]s, and its [`PackEngine::run`] **is** the lifecycle —
/// the only one in the crate. Each preset entry is one line
/// (`Pipeline::single(self).run(targets, max_loops)`), so a preset run and a
/// one-stage pipeline are the same code, not two implementations to keep in
/// step.
///
/// # The five parts of a run
///
/// 1. **Validate.** No targets / an empty molecule / each factory's
///    [`StageFactory::validate_targets`] — named, before anything is built.
/// 2. **Resolve space and build state.** Global restraints are broadcast onto
///    every target, the density / box / cell declaration resolves into one
///    packing space, and the [`PackContext`](crate::PackContext) is built and
///    wrapped — with a zeroed rigid placement vector — into one
///    [`PackState`].
/// 3. **Check the chain.** Every factory's stages are resolved, then three
///    refusals fire *before any handler is notified and before any stage
///    runs*: an empty chain ([`PackError::NoStages`]), a factory carrying
///    non-default shared settings ([`PackError::PresetSettingsInsidePipeline`]),
///    and a stage whose [`Requires`](crate::Requires) cannot hold where it
///    sits ([`PackError::StageOrder`]).
/// 4. **Run each stage.** `on_start` opens the handler bracket once for the
///    whole run; then, per stage: the geometry cache is invalidated,
///    `on_stage_start` fires, the stage runs, the placement marker advances
///    by the stage's [`Guarantees`](crate::Guarantees) (declared, never
///    inspected), `on_stage_end` fires, `degraded` accumulates. A handler
///    asking to stop ends the run there: later stages never start, and the
///    verdict is honestly `converged == false`. A stage that *fails* returns
///    a named [`PackError`] instead, which propagates on the same path as
///    part 1's refusals: no `on_stage_end` for it, no `on_finish` for the
///    run, no half-built result.
/// 5. **Assemble.** The lab-frame coordinates are rebuilt from the rigid view
///    with [`RigidView::write_xcart`](crate::RigidView::write_xcart),
///    `on_finish` closes the bracket, and the frame plus the placement
///    solution become the [`State`].
///
/// # Handlers are adopted; settings are refused
///
/// A preset handed to [`Pipeline::with_stage`] may carry two things. Its
/// handlers are **adopted** — appended to the pipeline's set in stage order,
/// observing the whole run, because observing across stages is exactly a
/// handler's semantics. Its shared [`PackSettings`] are **refused by name**:
/// tolerance, precision, seed and the cell are one ruler, and two stages
/// each holding one would leave the objective with no single ruler.
/// [`Pipeline::single`] is the other case — it *adopts* the engine's
/// settings, which is what makes a preset's own `run` legal, and why a
/// factory reads shared knobs from [`EngineSetup::settings`], never its own.
///
/// # Where the stage index comes from
///
/// A [`Stage`] builds its `StepInfo` with `StageInfo { index: 0, total: 1 }`
/// — correct when it is the whole run, and all it can know otherwise. The
/// **position** is the pipeline's fact, so the pipeline owns it: every
/// handler is wrapped once by `bracket.rs`'s private `StageTagger`, which
/// overwrites `StepInfo.stage` from a shared slot updated before each stage;
/// neither [`Stage`] nor [`PackState`] learns where it sits.
///
/// # Combinators: what they cost, and how to take them back out
///
/// [`with_repeat`](Pipeline::with_repeat) and
/// [`with_guarded`](Pipeline::with_guarded) build the two stages
/// the pipeline's combinators define, and both go through
/// [`with_stage`](Pipeline::with_stage), so handler adoption and the
/// settings refusal have one spelling. Their first real consumer is the
/// `dg-refine` recipe: "connect ⇄ refine" alternates to convergence
/// (`Repeat { Until::Converged }`), and its ring closure needs a guarded
/// retry (`Guarded { RingClosed, Rerun { max } }`). They rest on two earlier
/// promises — the re-entrancy contract on [`Stage::run`]
/// (a repeated stage keeps its own configuration) and the `Placed::All`
/// continuation below (pass *n+1* continues pass *n* instead of re-running
/// `initial()`). If that consumer is ever dropped, the whole feature comes
/// out as a unit: `combinators.rs`, `src/invariant.rs`, these two builders
/// and `PackError::InvariantViolated` — nothing else depends on them.
///
/// # The cache boundary
///
/// [`PackState::invalidate_geometry_cache`] is called before **every** stage.
/// A pipeline reuses one context, so the previous stage's last evaluation
/// leaves the geometry cache *hot*, while the hand-written `with_restart`
/// spelling starts from a fresh context and therefore a *cold* one. A hit
/// changes the summation path in `objective.rs` (it skips the cell reset and
/// the molecule expansion, accumulating the constraint values from `xcart`
/// instead), so `frest` would no longer agree bit for bit between the two
/// spellings. Invalidating at each boundary starts both cold; before the
/// first stage the context is new and the call is a no-op.
///
/// # One verdict, one `Placements`
///
/// `fdist` / `frest` are read off the **state** after the last stage ran —
/// never from a `StageOutcome`, which carries no verdict by design — and
/// `converged` is the last stage's own flag AND those two numbers under the
/// run's precision. The placement solution is likewise taken from the state's
/// rigid-view slot with no branch on which algorithm produced it, which rests
/// on the seam's writeback contract: **every stage leaves the state's rigid
/// view valid** (GENCAN writes its `x`; both growth drivers capture from
/// `xcart`). The bitwise continuity `with_restart` promises therefore holds
/// for a run whose *last* stage maintained that view — which every stage in
/// this crate does.
pub struct Pipeline {
    settings: PackSettings,
    handlers: Vec<Box<dyn Handler>>,
    factories: Vec<Box<dyn StageFactory>>,
}

impl Default for Pipeline {
    fn default() -> Self {
        Self::new()
    }
}

impl Pipeline {
    /// An empty pipeline: default shared settings, no handlers, no stages.
    /// Running it is [`PackError::NoStages`] — a named error, never a run
    /// that quietly hands the input back.
    pub fn new() -> Self {
        Self {
            settings: PackSettings::default(),
            handlers: Vec::new(),
            factories: Vec::new(),
        }
    }

    /// One engine as a whole run — the spelling every preset's `run` uses.
    ///
    /// Unlike [`with_stage`](Self::with_stage) this **adopts** the engine's
    /// shared settings — it is the only stage source, so its ruler is the
    /// run's — as well as its handlers. The engine is left holding the
    /// defaults, so the "no second ruler" check reads the same on both
    /// spellings and the run's knobs reach it through
    /// [`EngineSetup::settings`].
    pub fn single(mut engine: impl PackEngine + 'static) -> Self {
        let settings = std::mem::take(engine.settings_mut());
        let handlers = engine.take_handlers();
        Self {
            settings,
            handlers,
            factories: vec![Box::new(engine)],
        }
    }

    /// Append a stage source to the chain.
    ///
    /// The stage's handlers are **adopted** here, in stage order, and go on
    /// to observe the whole run. Its shared [`PackSettings`] are refused when
    /// they are not the defaults: the refusal
    /// ([`PackError::PresetSettingsInsidePipeline`], naming the knob) surfaces
    /// from [`run`](PackEngine::run) — this builder returns `Self` and has
    /// nowhere to put a `Result`. Declare shared knobs on the pipeline.
    pub fn with_stage(mut self, mut stage: impl StageFactory + 'static) -> Self {
        self.handlers.extend(stage.take_handlers());
        self.factories.push(Box::new(stage));
        self
    }

    /// Append `body` as one stage that runs it until `until` is met.
    ///
    /// Semantics — `Passes(0)`, the unbounded `Converged`, what a pass owes
    /// its stages — are on [`Until`] and in the `combinators` module docs.
    pub fn with_repeat(self, body: Vec<Box<dyn StageFactory>>, until: Until) -> Self {
        self.with_stage(RepeatFactory { body, until })
    }

    /// Append `stage` as one stage whose exit must satisfy `invariants`,
    /// answering a violation by `on_violation` — a rerun of the same stage
    /// or a named failure, never another algorithm (see [`OnViolation`]).
    pub fn with_guarded(
        self,
        stage: impl StageFactory + 'static,
        invariants: Vec<Box<dyn Invariant>>,
        on_violation: OnViolation,
    ) -> Self {
        self.with_stage(GuardedFactory {
            inner: Box::new(stage),
            invariants,
            on_violation,
        })
    }
}

impl StageFactory for Pipeline {
    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        for f in &self.factories {
            f.validate_targets(targets)?;
        }
        Ok(())
    }
    fn settings(&self) -> &PackSettings {
        &self.settings
    }
    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        std::mem::take(&mut self.handlers)
    }

    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        let mut out = Vec::new();
        for f in self.factories.iter_mut() {
            out.extend(f.stages(setup)?);
        }
        Ok(out)
    }
}

/// [`PackEngine::run`]'s five parts, split out (module doc numbering).
impl Pipeline {
    fn validate(&self, targets: &[Target]) -> Result<(), PackError> {
        if targets.is_empty() {
            return Err(PackError::NoTargets);
        }
        if let Some((i, _)) = targets.iter().enumerate().find(|(_, t)| t.natoms() == 0) {
            return Err(PackError::EmptyMolecule(i));
        }
        crate::assemble::check_templates(targets)?;
        for f in &self.factories {
            f.validate_targets(targets)?;
        }
        Ok(())
    }
    fn resolve_stages(
        factories: &mut [Box<dyn StageFactory>],
        setup: &EngineSetup<'_>,
    ) -> Result<Vec<Box<dyn Stage>>, PackError> {
        // Per factory, not flattened yet: the settings refusal has to name
        // the stage the offending factory produces.
        let mut per_factory: Vec<Vec<Box<dyn Stage>>> = Vec::new();
        for f in factories.iter_mut() {
            per_factory.push(f.stages(setup)?);
        }
        if per_factory.iter().all(|s| s.is_empty()) {
            return Err(PackError::NoStages);
        }
        for (f, produced) in factories.iter().zip(per_factory.iter()) {
            if let Some(knob) = f.settings().first_non_default_knob() {
                return Err(PackError::PresetSettingsInsidePipeline {
                    stage: produced.first().map_or("pipeline", |s| s.name()),
                    knob,
                });
            }
        }
        let stages: Vec<Box<dyn Stage>> = per_factory.into_iter().flatten().collect();

        let mut placed = Placed::None;
        for s in stages.iter() {
            if s.requires().placed == Placed::All && placed == Placed::None {
                return Err(PackError::StageOrder {
                    stage: s.name(),
                    needs: "placed: all",
                });
            }
            placed = s.guarantees().placed;
        }
        Ok(stages)
    }
    /// Open the handler bracket, then run every stage in order, keeping the
    /// shared stage position current before each one so `StageTagger`
    /// (`bracket.rs`) stamps the right identity on the `on_step` events a
    /// stage emits. Per stage: invalidate the geometry cache, fire
    /// `on_stage_start`/`on_stage_end`, `set_placed` by the stage's declared
    /// guarantee, and stop early the moment a handler asks for one. Returns
    /// `(last_converged, degraded)`, or the error a stage failed with —
    /// which skips its `on_stage_end` and the run's `on_finish`. This loop,
    /// not the bracket around it, is the heart of the lifecycle.
    fn run_stages(
        state: &mut PackState,
        stages: &mut [Box<dyn Stage>],
        setup: &EngineSetup<'_>,
        space: &ResolvedSpace,
        budget: &Budget,
        own_handlers: Vec<Box<dyn Handler>>,
        handlers: &mut Vec<Box<dyn Handler>>,
    ) -> Result<(bool, usize), PackError> {
        let (tagged, position) = open_bracket(own_handlers, stages, setup, space, budget);
        *handlers = tagged;
        state.ctx_mut().ntotmol = setup.ntotmol_free;

        let total = stages.len();
        let mut degraded = 0usize;
        let mut last_converged = false;
        for (index, stage) in stages.iter_mut().enumerate() {
            state.invalidate_geometry_cache();
            let info = StageInfo {
                index,
                total,
                name: stage.name(),
            };
            for h in handlers.iter_mut() {
                h.on_stage_start(&info);
            }
            *position.lock().expect("stage position mutex") = info;

            let outcome = stage.run(state, setup.targets, budget, handlers)?;
            // By what the stage declared, never by inspecting the result.
            state.set_placed(stage.guarantees().placed);

            for h in handlers.iter_mut() {
                h.on_stage_end(&info, &outcome, state.ctx());
            }
            degraded += outcome.degraded;
            last_converged = outcome.converged;
            if handlers.iter().any(|h| h.should_stop()) {
                last_converged = false;
                break;
            }
        }
        Ok((last_converged, degraded))
    }
    fn assemble(
        mut state: PackState,
        handlers: &mut [Box<dyn Handler>],
        setup: &EngineSetup<'_>,
        space: &ResolvedSpace,
        outcome: (bool, usize),
        precision: F,
    ) -> Result<State, PackError> {
        let (last_converged, degraded) = outcome;
        {
            let ctx = state.ctx_mut();
            for itype in 0..setup.ntype_with_fixed {
                ctx.comptype[itype] = true;
            }
            ctx.ntotmol = setup.ntotmol_free;
        }
        let (mut sys, view) = state.into_parts();
        view.write_xcart(&mut sys);
        close_bracket(handlers, &sys);
        // The verdict is read off the state the last stage left behind — the
        // shared objective's own numbers, never a stage's self-report.
        let (fdist, frest) = (sys.fdist, sys.frest);
        let converged = last_converged && fdist < precision && frest < precision;

        // The placement solution, verbatim, for cross-entry seeding: the frame
        // below is derived VIEW data — re-deriving (coor, rigid) from it would
        // recompute COMs and break bitwise continuity. Every stage leaves this
        // slot valid, so there is no branch on which one wrote it.
        let placements = result::Placements {
            rigid: view.clone(),
            coor: sys.coor[..setup.ntotat_free].to_vec(),
            copy_atoms: setup
                .targets
                .iter()
                .filter(|t| t.fixed_at.is_none())
                .flat_map(|t| std::iter::repeat_n(t.natoms(), t.count))
                .collect(),
            cell: sys.simbox.clone(),
        };

        let xcart = std::mem::take(&mut sys.xcart);
        let positions = positions_in_target_order(setup.targets, &xcart, setup.ntotat_free);
        let intra = IntraResidual::from_targets(setup.targets, &positions, &sys.simbox);
        let mut frame = crate::assemble::assemble_frame(setup.targets, &positions)?;
        // The resolved cell — a periodic box, a declared lattice or a seed's
        // inherited cell — is user-stated geometry and belongs on the output
        // frame; the fallback box drawn around the atoms (`sys.simbox` when no
        // cell was declared) is not.
        frame.simbox = space.cell.clone();

        Ok(State {
            frame,
            placements,
            fdist,
            intra,
            frest,
            converged,
            degraded,
        })
    }
}

impl PackEngine for Pipeline {
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>> {
        &mut self.handlers
    }

    fn run(mut self, targets: &[Target], max_loops: usize) -> Result<State, PackError> {
        // ── 1. Validate ───────────────────────────────────────────────────
        self.validate(targets)?;

        // ── 2. Space, context, state ──────────────────────────────────────
        let broadcast = broadcast_global_restraints(targets, &self.settings.global_restraints);
        let targets: &[Target] = &broadcast;
        let space = resolve_pack_space(
            self.settings.density,
            self.settings.periodic_box,
            self.settings.cell.clone(),
            targets,
        )?;

        let knobs = ContextKnobs {
            tolerance: self.settings.tolerance(),
            short_tolerance: self.settings.short_tolerance,
            parallel_eval: self.settings.parallel_eval,
        };
        let built = build_context(&knobs, targets)?;
        let mut state = PackState::new(built.sys, built.ntotmol_free);

        let setup = EngineSetup {
            settings: &self.settings,
            targets,
            cell: space.cell.clone(),
            maxmove_per_type: &built.maxmove_per_type,
            ntype: built.ntype,
            ntype_with_fixed: built.ntype_with_fixed,
            ntotmol_free: built.ntotmol_free,
            ntotat: built.ntotat,
            ntotat_free: built.ntotat_free,
        };

        // ── 3. Check the chain (before any handler is notified) ───────────
        let mut stages = Self::resolve_stages(&mut self.factories, &setup)?;

        // ── 4. Run the stages, bracketed by the handlers ──────────────────
        let precision = self.settings.precision();
        let budget = Budget::new(max_loops, precision);
        let own_handlers = std::mem::take(&mut self.handlers);
        let mut handlers: Vec<Box<dyn Handler>> = Vec::new();
        let outcome = Self::run_stages(
            &mut state,
            &mut stages,
            &setup,
            &space,
            &budget,
            own_handlers,
            &mut handlers,
        )?;

        // ── 5. Rebuild, close the bracket, assemble ───────────────────────
        Self::assemble(state, &mut handlers, &setup, &space, outcome, precision)
    }
}
