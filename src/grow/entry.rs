//! `CbmcGrow` — the continuum configurational-bias chain-growth entry.

use crate::Handler;
use crate::PackError;
use crate::Stage;
use crate::Target;
use crate::entry::{PackSettings, State};
use crate::grow::GrowStage;
use crate::grow::config::GrowConfig;
use crate::grow::prior::{AnglePrior, TorsionPrior};
use crate::grow::validate_grow_cell;
use crate::pipeline::{EngineSetup, PackEngine, Pipeline, StageFactory};
use molrs::op::F;

/// Configurational-bias chain growth (CBMC-style constructive packing) as
/// its own entry.
///
/// The torsion prior is mandatory — it decides the grown chains' statistics
/// (spec Domain basis §5.1) — so it is the one constructor argument.
/// Shared knobs come from [`PackEngine`]; the growth knobs forward to the
/// underlying [`GrowConfig`].
///
/// `max_loops` — the second argument to [`PackEngine::run`] — is read as an
/// allowance of *passes over a chain*, not of individual moves. Growing a chain
/// costs one round per stage, so the driver's round loop is capped at
/// `max_loops × (the longest chain's number of steps + 1)` rounds. Hitting that
/// cap force-completes the unfinished chains instead of spinning forever, which
/// grows `degraded` and leaves `converged == false`.
///
/// The entry reports its outcome honestly: on non-convergence `degraded`
/// and `converged` say so and nothing else runs — no hidden second
/// algorithm (engine-entry-split 门槛 2). For the rigid push-off, chain
/// explicitly: feed the same free targets to
/// [`GenCanPack::with_restart`](crate::GenCanPack::with_restart) with this
/// run's result (placement-seeding spec).
pub struct CbmcGrow {
    settings: PackSettings,
    handlers: Vec<Box<dyn Handler>>,
    config: GrowConfig,
}

impl CbmcGrow {
    pub fn new(torsion_prior: TorsionPrior) -> Self {
        Self::from_config(GrowConfig::new(torsion_prior))
    }

    /// Build the entry around an existing [`GrowConfig`].
    pub fn from_config(config: GrowConfig) -> Self {
        Self {
            settings: PackSettings::default(),
            handlers: Vec::new(),
            config,
        }
    }

    /// Torsion candidates proposed per growth step (CBMC `k`).
    pub fn with_trials(mut self, trials: usize) -> Self {
        self.config = self.config.with_trials(trials);
        self
    }
    /// Rosenbluth inverse temperature applied to the crowding penalty.
    pub fn with_selectivity(mut self, beta: F) -> Self {
        self.config = self.config.with_selectivity(beta);
        self
    }
    /// Soft-shell width in Å beyond the hard core.
    pub fn with_soft_shell(mut self, width: F) -> Self {
        self.config = self.config.with_soft_shell(width);
        self
    }
    /// Steps retracted on a dead end.
    pub fn with_retract(mut self, steps: usize) -> Self {
        self.config = self.config.with_retract(steps);
        self
    }
    /// Tail regrowth cadence and window (`every = 0` disables).
    pub fn with_relax(mut self, every: usize, window: usize) -> Self {
        self.config = self.config.with_relax(every, window);
        self
    }
    /// Cumulative dead ends on a chain before that chain's hard core softens by one rung
    /// (one rung multiplies the dimensionless hard-core scale by
    /// [`GrowConfig::SOFTEN_RUNG`](crate::grow::config::GrowConfig::SOFTEN_RUNG)).
    ///
    /// Forwards to [`GrowConfig::with_soften_after`]. The clock is **cumulative**
    /// dead ends on that chain (`deadends_total` versus the `rungs_earned`
    /// watermark), **not** consecutive. A successful commit does not reset
    /// the softening counter.
    ///
    /// `rung_due` consumes this cadence; `force_due` consumes
    /// `2 × soften_after` consecutive streak at the floor; `retract_depth`
    /// reads the streak and the separate retract knob, never this counter.
    /// The Python wheel and `docs/python/api-reference.md` still describe this
    /// knob as consecutive; that page is left stale on purpose until
    /// special-bonds-06.
    pub fn with_soften_after(mut self, attempts: usize) -> Self {
        self.config = self.config.with_soften_after(attempts);
        self
    }
    /// Softening floor for the dimensionless hard-core scale.
    ///
    /// Forwards to [`GrowConfig::with_min_hard_scale`]. `1.0` is full declared
    /// contact (`radius_i + radius_j`); the ladder walks the scale down by
    /// [`GrowConfig::SOFTEN_RUNG`](crate::grow::config::GrowConfig::SOFTEN_RUNG)
    /// per rung to this floor (default 0.8, Auhl's 0.8σ floor, where σ
    /// is the excluded-volume / bead diameter).
    pub fn with_min_hard_scale(mut self, scale: F) -> Self {
        self.config = self.config.with_min_hard_scale(scale);
        self
    }
    /// Placement-angle prior (all-atom template default; WLC for CG).
    pub fn with_angle_prior(mut self, prior: AnglePrior) -> Self {
        self.config = self.config.with_angle_prior(prior);
        self
    }
    /// Grow one chain to completion before starting the next.
    pub fn with_serial(mut self, serial: bool) -> Self {
        self.config = self.config.with_serial(serial);
        self
    }
    /// Seed chains in empty field cells (cavity seeding).
    pub fn with_void_bias(mut self, void_bias: bool) -> Self {
        self.config = self.config.with_void_bias(void_bias);
        self
    }
}

impl StageFactory for CbmcGrow {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        for (i, t) in targets.iter().enumerate() {
            crate::grow::validate_template(t.template.as_ref())
                .map_err(|source| PackError::Grow { target: i, source })?;
            if t.fixed_at.is_some() {
                return Err(PackError::Grow {
                    target: i,
                    source: crate::grow::GrowError::FixedTarget,
                });
            }
        }
        Ok(())
    }

    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        std::mem::take(self.handlers_mut())
    }

    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        // Growth needs a box to grow into — a named error, resolved here so
        // it is reported before any stage runs. The stage installs it (and
        // the matching cell grid) at the top of its own run.
        let stage = GrowStage::from_targets(setup.targets, &self.config, setup.settings.seed())
            .map_err(|(target, source)| PackError::Grow { target, source })?;
        let cell = validate_grow_cell(setup.cell.clone())
            .map_err(|source| PackError::Grow { target: 0, source })?;
        Ok(vec![Box::new(
            stage.with_resolved_cell(cell, setup.settings.discale()),
        )])
    }
}

impl PackEngine for CbmcGrow {
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>> {
        &mut self.handlers
    }

    fn run(self, targets: &[Target], max_loops: usize) -> Result<State, PackError> {
        Pipeline::single(self).run(targets, max_loops)
    }
}
