//! [`PackSettings`] — the knobs every packing run shares.
//!
//! The knobs the shared infrastructure reads (contact tolerance, precision,
//! seed, the box or cell declaration, global restraints, screen logging via
//! [`LogSpec`]). One ruler per run: algorithm-specific knobs live on their
//! own engine type, never here. Two siblings complete what every run needs
//! and no algorithm owns: [`pack_space`](crate::pack_space) resolves a
//! density / periodic box / cell declaration into the one space the run
//! packs into and broadcasts global restraints onto every target, and
//! [`state`](crate::state) holds [`State`](crate::State) and the verbatim
//! placement solution it carries, which is what makes one run continuable
//! from another.
//!
//! What is deliberately *not* here: the engines themselves — [`GencanPack`](crate::GencanPack)
//! lives with the rigid-body family, [`CbmcGrow`](crate::CbmcGrow) and
//! [`LatticeGrow`](crate::LatticeGrow) with
//! growth — and the lifecycle that drives them, which
//! belongs to the type that owns it, [`Pipeline`](crate::Pipeline). This
//! module names neither: settings and space are read by the lifecycle, they
//! do not run it. The dependency arrow points one way only.

use molrs::op::F;

use crate::LogLevel;
use crate::pack_space::{CellDecl, PeriodicSpec};

/// Built-in screen logging: detail level + print cadence.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LogSpec {
    pub(crate) level: LogLevel,
    pub(crate) frequency: usize,
}

impl Default for LogSpec {
    fn default() -> Self {
        Self {
            level: LogLevel::Quiet,
            frequency: 1,
        }
    }
}

/// The knobs every engine shares — consumed by the lifecycle and the shared
/// infrastructure (space resolution, system construction), never by one
/// algorithm alone. Algorithm-specific knobs live on their engine.
#[derive(Debug, Clone, Default)]
pub struct PackSettings {
    pub(crate) tolerance: Option<F>,
    pub(crate) precision: Option<F>,
    pub(crate) discale: Option<F>,
    pub(crate) seed: Option<u64>,
    pub(crate) parallel_eval: bool,
    pub(crate) short_tolerance: Option<(F, F)>,
    pub(crate) periodic_box: Option<PeriodicSpec>,
    pub(crate) density: Option<F>,
    pub(crate) cell: Option<CellDecl>,
    pub(crate) log: LogSpec,
    pub(crate) global_restraints: Vec<std::sync::Arc<dyn crate::AtomRestraint>>,
}

impl PackSettings {
    /// Resolved contact tolerance (default 2.0 Å).
    pub fn tolerance(&self) -> F {
        self.tolerance.unwrap_or(2.0)
    }
    /// Resolved convergence precision (default 0.01).
    pub fn precision(&self) -> F {
        self.precision.unwrap_or(0.01)
    }
    /// Resolved initial radius up-scaling (default 1.1).
    pub fn discale(&self) -> F {
        self.discale.unwrap_or(1.1)
    }
    /// Resolved RNG seed (default 1_234_567).
    pub fn seed(&self) -> u64 {
        self.seed.unwrap_or(1_234_567)
    }

    /// The name of the first knob this set holds that is not the default, or
    /// `None` when every knob is untouched.
    ///
    /// The question a multi-stage run has to ask of a preset it was handed:
    /// these settings are *shared*, so a stage that carries its own would
    /// give the run a second ruler. The answer is a field name, which is what
    /// lets the refusal tell the caller exactly what to move.
    ///
    /// The body destructures [`PackSettings`] **completely** — no `..` — so a
    /// knob added later fails to compile here instead of slipping silently
    /// past the check.
    pub fn first_non_default_knob(&self) -> Option<&'static str> {
        let default = PackSettings::default();
        let PackSettings {
            tolerance,
            precision,
            discale,
            seed,
            parallel_eval,
            short_tolerance,
            periodic_box,
            density,
            cell,
            log,
            global_restraints,
        } = self;

        if *tolerance != default.tolerance {
            return Some("tolerance");
        }
        if *precision != default.precision {
            return Some("precision");
        }
        if *discale != default.discale {
            return Some("discale");
        }
        if *seed != default.seed {
            return Some("seed");
        }
        if *parallel_eval != default.parallel_eval {
            return Some("parallel_eval");
        }
        if *short_tolerance != default.short_tolerance {
            return Some("short_tolerance");
        }
        if *periodic_box != default.periodic_box {
            return Some("periodic_box");
        }
        if *density != default.density {
            return Some("density");
        }
        // `CellDecl` carries floats a caller cannot be expected to reproduce
        // exactly, and the default is "no declaration", so presence is the
        // whole question.
        if cell.is_some() {
            return Some("cell");
        }
        if log.level != default.log.level || log.frequency != default.log.frequency {
            return Some("log");
        }
        // Restraint objects are `dyn` and not comparable; the default is the
        // empty set, so emptiness is the whole question.
        if !global_restraints.is_empty() {
            return Some("global_restraints");
        }
        None
    }
}
