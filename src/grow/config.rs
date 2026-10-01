//! Leaf configuration for growth.
//!
//! The `CbmcGrow` entry carries a [`GrowConfig`], so this file must stay a
//! **leaf**: importing `target` / `entry` / `context` from here would close
//! a dependency cycle between the entry layer and the growth module. Only
//! sibling leaves (`prior`) and `molrs` types are allowed.

use molrs::BondDistanceWeights;
use molrs::types::F;

use crate::grow::prior::{AnglePrior, TorsionPrior};

/// Configuration for the chain-growth solver.
///
/// The torsion prior has **no default** — it decides the grown chains'
/// statistics (spec Domain basis §5.1), so the caller must state it, even if
/// the statement is `TorsionPrior::Uniform`. Everything else defaults.
/// The intramolecular skip table lives on `Target.special_bonds`.
#[derive(Debug, Clone)]
pub struct GrowConfig {
    pub(crate) torsion_prior: TorsionPrior,
    /// Torsion candidates proposed per growth step (CBMC `k`).
    pub(crate) trials: usize,
    /// Rosenbluth inverse temperature `β` applied to the crowding penalty.
    pub(crate) selectivity: F,
    /// Soft-shell width in Å beyond the hard core: allowed but charged.
    pub(crate) soft_shell: F,
    /// Steps retracted on a dead end (recoil at feeler depth 1).
    pub(crate) retract: usize,
    /// Regrow each chain's tail every this many rounds (`0` disables).
    pub(crate) relax_every: usize,
    /// Tail length (in steps) of the periodic regrowth.
    pub(crate) relax_window: usize,
    /// Cumulative dead ends on a chain per softening rung (`rung_due`
    /// watermark; one rung multiplies the hard-core scale by
    /// [`Self::SOFTEN_RUNG`]). Not consecutive: commits do not reset this
    /// cadence.
    pub(crate) soften_after: usize,
    /// Floor for the dimensionless hard-core scale. `1.0` is full declared
    /// contact (`radius_i + radius_j`); the softening ladder walks the scale
    /// down by [`Self::SOFTEN_RUNG`] per rung, never below this floor. Default
    /// 0.8 — Auhl's
    /// 0.8σ push-off bound, where σ is the excluded-volume (bead) diameter
    /// (Auhl et al. 2003). Contacts tighter than that floor are a poor
    /// starting point for a subsequent excluded-volume push-off.
    pub(crate) min_hard_scale: F,
    /// Placement-angle prior (`Template` = all-atom default; `Wlc` for CG).
    pub(crate) angle_prior: AnglePrior,
    /// Serial scheduling: grow one chain to completion before starting the
    /// next, instead of advancing every pending chain per round. Each new
    /// chain then threads a *finished* matrix rather than a box of partial
    /// chains (Packmol-style sequential insertion, but constructive).
    pub(crate) serial: bool,
    /// Void-biased seeding: draw seed (and forced-placement) anchors from
    /// empty field cells instead of uniformly over the box — the
    /// constructive half of Mezei's cavity-biased insertion (Mol. Phys. 40,
    /// 901 (1980)); no detailed-balance correction is owed because the
    /// solver constructs, it does not sample an ensemble.
    pub(crate) void_bias: bool,
}

impl GrowConfig {
    /// One softening rung. `hard_scale` is multiplied by this and never
    /// falls below [`Self::min_hard_scale`].
    pub const SOFTEN_RUNG: F = 0.97;

    /// A growth configuration with the mandatory torsion prior and defaults
    /// for everything else.
    pub fn new(torsion_prior: TorsionPrior) -> Self {
        Self {
            torsion_prior,
            trials: 12,
            selectivity: 2.0,
            soft_shell: 1.0,
            retract: 10,
            relax_every: 25,
            relax_window: 6,
            soften_after: 50,
            min_hard_scale: 0.8,
            angle_prior: AnglePrior::Template,
            serial: false,
            void_bias: false,
        }
    }

    /// Grow chains one at a time (each chain completes before the next
    /// starts) instead of the default round-robin over all pending chains.
    pub fn with_serial(mut self, serial: bool) -> Self {
        self.serial = serial;
        self
    }

    /// Draw seed anchors from empty field cells (cavity seeding) instead of
    /// uniformly over the box. The probe stays the arbiter — this only
    /// steers trials toward voids.
    pub fn with_void_bias(mut self, void_bias: bool) -> Self {
        self.void_bias = void_bias;
        self
    }

    /// The torsion prior this configuration was built with.
    ///
    /// **Rust-only:** not mirrored in the Python wheel, whose `GrowConfig`
    /// is write-only (constructor + builders).
    pub fn torsion_prior(&self) -> &TorsionPrior {
        &self.torsion_prior
    }

    /// Torsion candidates proposed per growth step.
    pub fn with_trials(mut self, trials: usize) -> Self {
        self.trials = trials.max(1);
        self
    }

    /// Rosenbluth inverse temperature applied to the crowding penalty.
    pub fn with_selectivity(mut self, beta: F) -> Self {
        self.selectivity = beta.max(0.0);
        self
    }

    /// Soft-shell width in Å beyond the hard core.
    pub fn with_soft_shell(mut self, width: F) -> Self {
        self.soft_shell = width.max(0.0);
        self
    }

    /// Steps retracted on a dead end.
    pub fn with_retract(mut self, steps: usize) -> Self {
        self.retract = steps.max(1);
        self
    }

    /// Regrow each chain's last `window` steps every `every` rounds
    /// (`every = 0` disables the periodic regrowth).
    pub fn with_relax(mut self, every: usize, window: usize) -> Self {
        self.relax_every = every;
        self.relax_window = window.max(1);
        self
    }

    /// Cumulative dead ends on a chain before the hard core softens by one rung
    /// (one rung multiplies the dimensionless hard-core scale by [`Self::SOFTEN_RUNG`]).
    ///
    /// The clock is **cumulative** dead ends on that chain (`deadends_total`
    /// versus the `rungs_earned` watermark), **not** consecutive. A successful
    /// commit does not reset the softening counter. Clamped ≥ 1; default 50.
    ///
    /// The growth driver keeps three private readers. `rung_due` consumes this
    /// cadence (`deadends_total` against `rungs_earned`) to earn a softening
    /// rung. `force_due` consumes `2 × soften_after` consecutive
    /// `deadend_streak` at the `min_hard_scale` floor to force-place a wedged
    /// chain. `retract_depth` reads the consecutive streak and the separate
    /// [`with_retract`](Self::with_retract) knob, never this counter.
    ///
    /// The Python wheel and `docs/python/api-reference.md` still describe this
    /// knob as consecutive; that page is left stale on purpose until
    /// special-bonds-06.
    pub fn with_soften_after(mut self, attempts: usize) -> Self {
        self.soften_after = attempts.max(1);
        self
    }

    /// Softening floor for the dimensionless hard-core scale.
    ///
    /// `hard_scale` is dimensionless: `1.0` is full declared contact
    /// (`radius_i + radius_j`); the ladder walks it down by [`Self::SOFTEN_RUNG`]
    /// per rung to this floor. Clamped to `[0.0, 1.0]`; default 0.8 is Auhl's 0.8σ
    /// push-off bound, where σ is the excluded-volume (bead) diameter
    /// (Auhl et al. 2003). Contacts tighter than that floor are a poor
    /// starting point for a subsequent excluded-volume push-off.
    pub fn with_min_hard_scale(mut self, scale: F) -> Self {
        self.min_hard_scale = scale.clamp(0.0, 1.0);
        self
    }

    /// Placement-angle prior. Defaults to [`AnglePrior::Template`] (angles
    /// copied verbatim — the all-atom behavior); CG templates use
    /// [`AnglePrior::Wlc`] to control persistence (spec §5.5).
    pub fn with_angle_prior(mut self, prior: AnglePrior) -> Self {
        self.angle_prior = prior;
        self
    }
}

/// First table slot whose weight is neither `0.0` nor `1.0`.
///
/// Slot 0 is the 1-2 weight. `None` means the table is binary and growth
/// may compile it as a skip set.
pub(crate) fn binary_violation(weights: &BondDistanceWeights) -> Option<(usize, F)> {
    weights
        .as_slice()
        .iter()
        .copied()
        .enumerate()
        .find(|&(_, w)| w != 0.0 && w != 1.0)
}

/// Rosenbluth crowding cap shared by propose and commit.
///
/// Both steps must use this value: a trial accepted against one shell and
/// scored against another is a different move than the one that was chosen.
pub(crate) fn crowding_cap(selectivity: F) -> F {
    60.0 / selectivity.max(0.1)
}
