//! Leaf configuration and error types for growth.
//!
//! The `CbmcGrow` entry carries a [`GrowConfig`], so this file must stay a
//! **leaf**: importing `target` / `entry` / `context` from here would close
//! a dependency cycle between the entry layer and the growth module. Only
//! sibling leaves (`prior`) and `molrs` types are allowed.

use std::fmt;

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
    /// watermark; one rung multiplies the hard-core scale by 0.97). Not
    /// consecutive: commits do not reset this cadence.
    pub(crate) soften_after: usize,
    /// Floor for the dimensionless hard-core scale. `1.0` is full declared
    /// contact (`radius_i + radius_j`); the softening ladder walks the scale
    /// down by 0.97 per rung, never below this floor. Default 0.8 — Auhl's
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
    /// (one rung multiplies the dimensionless hard-core scale by 0.97).
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
    /// (`radius_i + radius_j`); the ladder walks it down by 0.97 per rung
    /// to this floor. Clamped to `[0.0, 1.0]`; default 0.8 is Auhl's 0.8σ
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

/// Why a target cannot be grown.
///
/// Growth consumes the template's *chemistry* (its bond graph), so a target
/// that carries none is refused with a named error — never silently degraded
/// to rigid-body packing. The caller chooses the method per target.
///
/// Template-graph refusals stop at the first match, in this order:
/// `NoAtomsBlock → BondOutOfRange → NoBonds → TemplateTooSmall →
/// Disconnected → RingTemplate`.
#[derive(Debug, Clone)]
pub enum GrowError {
    /// The target was built without a template frame
    /// ([`Target::from_coords`][crate::Target::from_coords]), so there is no
    /// bond graph to grow from.
    MissingTemplate,
    /// The template has fewer than 3 atoms; growth needs a rigid seed of 3.
    TemplateTooSmall(usize),
    /// The template frame has no readable `atoms` block (`x` / `y` / `z` in Å).
    NoAtomsBlock,
    /// A bond names an endpoint outside the template.
    BondOutOfRange {
        /// First endpoint of the offending bond, as written in the frame.
        a: usize,
        /// Second endpoint of the offending bond, as written in the frame.
        b: usize,
        /// Atom count of the template (`xyz.len()`).
        n: usize,
    },
    /// The template frame carries no bonds (missing/empty graph).
    NoBonds,
    /// The template's bond graph does not connect all atoms.
    Disconnected,
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
    /// A special-bonds weight is neither 0 nor 1. Growth compiles a binary
    /// skip table; `index` is the 0-based slot (slot 0 is 1-2, slot 2 is 1-4).
    NonBinarySpecialBond {
        /// 0-based table slot.
        index: usize,
        /// The stored weight at `index`.
        weight: F,
    },
    /// The template's heavy-atom backbone does not fit the diamond lattice
    /// (the message names the offense). Tetrahedral heavy degree `1..=4` is
    /// accepted (linear is the `d = 2` degeneracy); degree `> 4`, a
    /// detached heavy, or a non-sp³ interior bond is named here.
    NonTetrahedralTemplate(String),
    /// The attached region contains no usable diamond site for this chain
    /// (empty Region ∩ lattice, or the tree cannot embed in the allowed
    /// subgraph).
    LatticeRegionEmpty,
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
            GrowError::NoAtomsBlock => {
                write!(f, "the template frame has no readable atoms block")
            }
            GrowError::BondOutOfRange { a, b, n } => write!(
                f,
                "bond ({a}, {b}) references an atom outside the template (natoms = {n})"
            ),
            GrowError::NoBonds => write!(
                f,
                "the template frame carries no bonds; growth needs the bond graph — pack \
                 this target with GenCanPack or supply connectivity"
            ),
            GrowError::Disconnected => {
                write!(f, "the template's bond graph does not connect all atoms")
            }
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
            GrowError::NonBinarySpecialBond { index, weight } => {
                let pair = index + 2;
                write!(
                    f,
                    "special-bonds 1-{pair} weight is {weight}; growth compiles a \
                     binary skip table (0 or 1). For all-atom explicit hydrogen, \
                     use Target::with_atom_radius rather than a fractional weight"
                )
            }
            GrowError::NonTetrahedralTemplate(msg) => write!(
                f,
                "the template's backbone does not fit the diamond lattice: {msg}"
            ),
            GrowError::FixedTarget => write!(
                f,
                "a fixed target cannot be grown: drop fixed_at or pack this \
                 target with GenCanPack"
            ),
            GrowError::LatticeRegionEmpty => write!(
                f,
                "LatticeGrow: the attached region contains no usable diamond \
                 site for this chain; enlarge the mesh or reduce the template"
            ),
        }
    }
}

impl std::error::Error for GrowError {}
