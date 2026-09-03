//! Leaf configuration and error types for growth.
//!
//! The `CbmcGrow` entry carries a [`GrowConfig`], so this file must stay a
//! **leaf**: importing `target` / `entry` / `context` from here would close
//! a dependency cycle between the entry layer and the growth module. Only
//! leaves (`prior`, and the crate-root [`topology`](crate::topology), which
//! itself imports nothing from this crate) and `molrs` types are allowed.

use std::fmt;

use molrs::types::F;

use crate::grow::prior::{AnglePrior, TorsionPrior};
use crate::topology::TopologyError;

/// Configuration for the chain-growth solver.
///
/// The torsion prior has **no default** — it decides the grown chains'
/// statistics (spec Domain basis §5.1), so the caller must state it, even if
/// the statement is `TorsionPrior::Uniform`. Everything else defaults.
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
    /// Consecutive dead ends at one step before the hard core is softened.
    pub(crate) soften_after: usize,
    /// Floor for the *hard-core scale*: the factor multiplying the pair
    /// contact distance (in Å) a candidate placement must clear. It starts at
    /// 1.0 — the full declared tolerance — and the softening ladder walks it
    /// down in 3 % rungs, never below this floor. Default 0.8, i.e. 80 % of
    /// tolerance: the push-off bound of Auhl et al. 2003, below which a
    /// downstream force-field relaxation can no longer pull contacts apart.
    pub(crate) min_hard_scale: F,
    /// Intramolecular exclusion depth in bonds (3 = 1-2/1-3/1-4, the
    /// all-atom convention; CG templates typically use 1 or 2).
    pub(crate) exclusion_depth: usize,
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
            exclusion_depth: 3,
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

    /// Consecutive dead ends at one step before the hard core softens.
    pub fn with_soften_after(mut self, attempts: usize) -> Self {
        self.soften_after = attempts.max(1);
        self
    }

    /// Softening floor for the hard-core scale — the factor multiplying the
    /// pair contact distance (in Å) a candidate placement must clear. 1.0 is
    /// the full declared tolerance, and the softening ladder walks the scale
    /// down in 3 % rungs but never past this floor. Clamped to `[0.0, 1.0]`;
    /// the default 0.8 (80 % of tolerance) is the push-off bound of Auhl et al.
    /// 2003, below which a downstream force-field relaxation can no longer pull
    /// contacts apart.
    pub fn with_min_hard_scale(mut self, scale: F) -> Self {
        self.min_hard_scale = scale.clamp(0.0, 1.0);
        self
    }

    /// Intramolecular exclusion depth in bonds.
    pub fn with_exclusion_depth(mut self, depth: usize) -> Self {
        self.exclusion_depth = depth;
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

/// Why a target cannot be grown.
///
/// Growth consumes the template's *chemistry* (its bond graph), so a target
/// that carries none is refused with a named error — never silently degraded
/// to rigid-body packing. The caller chooses the method per target.
#[derive(Debug, Clone)]
pub enum GrowError {
    /// The target was built without a template frame
    /// ([`Target::from_coords`][crate::Target::from_coords]), so there is no
    /// bond graph to grow from.
    MissingTemplate,
    /// The template has fewer than 3 atoms; growth needs a rigid seed of 3.
    ///
    /// The bond graph is read first, so a template that is both bondless and
    /// too small reports [`GrowError::Topology`] with `NoBonds`; the full order
    /// is `NoAtomsBlock → NoBonds → BondOutOfRange → TemplateTooSmall →
    /// Disconnected → RingTemplate`.
    TemplateTooSmall(usize),
    /// The template's bond graph could not be read, or does not qualify for
    /// growth. The bond graph is owned by [`crate::topology`]; growth wraps
    /// its error rather than restating the variants.
    ///
    /// Its `Display` is passed through verbatim, so the user sees the leaf's
    /// wording; match the variant to recover the [`TopologyError`] (there is no
    /// `Error::source` chain).
    Topology(TopologyError),
    /// No placed reference atom could be found while decomposing atom `.0`.
    NoReference(usize),
    /// Rotatable-bond perception failed.
    Perceive(String),
    /// A Grow target needs a box: neither a periodic box nor a cell (nor a
    /// density, once `with_density` lands) was declared.
    NoBox,
    /// The declared cell is not orthorhombic; the v1 overlap field only
    /// supports orthorhombic boxes (`triclinic-cell-downshift` lifts this).
    TriclinicCell,
    /// `fixed_at` combined with a growth entry — a fixed placement is by
    /// definition not grown.
    FixedTarget,
    /// The template's bond graph contains a cycle. Growth decomposes the
    /// template into a *tree* of internal coordinates; a ring bond would be
    /// silently dropped and the ring grown open — refused instead (lattice
    /// ring closure is its own future spec).
    RingTemplate,
    /// The template's heavy-atom backbone does not fit the diamond lattice:
    /// branched or non-tetrahedral (the message names the offense). Lattice
    /// growth v1 maps a linear sp³ backbone; branched trees are staged
    /// (lattice-growth-phase spec, 拓扑范围).
    NonTetrahedralTemplate(String),
}

impl fmt::Display for GrowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GrowError::MissingTemplate => write!(
                f,
                "the target has no template frame (built from bare coordinates); growth needs \
                 a bond graph — load the species from a file or build a frame with bonds, or \
                 pack this target with GenCanPack"
            ),
            GrowError::TemplateTooSmall(n) => write!(
                f,
                "the template has {n} atom(s); growth needs at least 3 — pack this \
                 target with GenCanPack"
            ),
            // Verbatim pass-through: the topology leaf owns this wording.
            GrowError::Topology(e) => write!(f, "{e}"),
            GrowError::NoReference(i) => write!(
                f,
                "no placed reference atom found while decomposing atom {i}"
            ),
            GrowError::Perceive(msg) => write!(f, "rotatable-bond perception failed: {msg}"),
            GrowError::NoBox => write!(
                f,
                "growth needs a box: declare with_periodic_box / with_cell (or a density) \
                 — the box is at its final volume from the first atom"
            ),
            GrowError::TriclinicCell => write!(
                f,
                "growth currently supports orthorhombic boxes only; declare an \
                 orthorhombic cell or periodic box"
            ),
            GrowError::RingTemplate => write!(
                f,
                "the template contains a ring: growth decomposes the bond graph \
                 into a tree, and a ring bond would be silently dropped; ring \
                 templates are refused"
            ),
            GrowError::NonTetrahedralTemplate(msg) => write!(
                f,
                "the template's backbone does not fit the diamond lattice: {msg}"
            ),
            GrowError::FixedTarget => write!(
                f,
                "a fixed target cannot be grown: drop fixed_at or pack this \
                 target with GenCanPack"
            ),
        }
    }
}

impl std::error::Error for GrowError {}
