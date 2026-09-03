//! `CbmcGrow` — the continuum configurational-bias chain-growth entry.

use crate::entry::{EngineSetup, PackEngine, PackSettings};
use crate::error::PackError;
use crate::grow::GrowthSolver;
use crate::grow::config::GrowConfig;
use crate::grow::prior::{AnglePrior, TorsionPrior};
use crate::grow::validate_grow_cell;
use crate::handler::Handler;
use crate::solver::Solver;
use crate::target::Target;
use molrs::types::F;

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
/// grows `softened` and leaves `converged == false`.
///
/// The entry reports its outcome honestly: on non-convergence `softened`
/// and `converged` say so and nothing else runs — no hidden second
/// algorithm (engine-entry-split 门槛 2). For the rigid push-off, chain
/// explicitly: feed the same free targets to
/// [`GenCanPack::seeded_from`](crate::GenCanPack::seeded_from) with this
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
    /// Consecutive dead ends before the hard core softens.
    pub fn with_soften_after(mut self, attempts: usize) -> Self {
        self.config = self.config.with_soften_after(attempts);
        self
    }
    /// Softening floor for the hard-core scale.
    pub fn with_min_hard_scale(mut self, scale: F) -> Self {
        self.config = self.config.with_min_hard_scale(scale);
        self
    }
    /// Intramolecular exclusion depth in bonds.
    pub fn with_exclusion_depth(mut self, depth: usize) -> Self {
        self.config = self.config.with_exclusion_depth(depth);
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

impl PackEngine for CbmcGrow {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }
    fn settings_mut(&mut self) -> &mut PackSettings {
        &mut self.settings
    }
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>> {
        &mut self.handlers
    }

    fn validate(&self, targets: &[Target]) -> Result<(), PackError> {
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

    fn prepare(
        &self,
        sys: &mut crate::context::PackContext,
        _x: &mut [F],
        setup: &EngineSetup<'_>,
    ) -> Result<(), PackError> {
        let simbox = validate_grow_cell(setup.cell.clone(), 0)?;
        let radmax = sys.radius.iter().cloned().fold(0.0 as F, F::max);
        crate::initial::install_simbox_and_grid(
            sys,
            simbox,
            radmax,
            self.settings.discale(),
            setup.ntotat_free,
        );
        Ok(())
    }

    fn solver(&mut self, setup: &EngineSetup<'_>) -> Result<Box<dyn Solver>, PackError> {
        let solver = GrowthSolver::from_targets(setup.targets, &self.config, self.settings.seed())
            .map_err(|(target, source)| PackError::Grow { target, source })?;
        Ok(Box::new(solver))
    }
}
