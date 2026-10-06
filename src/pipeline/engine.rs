//! What a run is made of: the resolved setup, the factory that turns it into
//! stages, and the engine surface a caller builds and runs.
//!
//! Two traits, one direction:
//!
//! * [`StageFactory`] — "given the resolved setup, what stages do you
//!   contribute?" Everything a [`Pipeline`](super::Pipeline) needs from the
//!   things it was handed: their target validation, their shared settings,
//!   the handlers they carry, and their stages. A preset entry implements it,
//!   and so does `Pipeline` itself, which is what makes pipelines nest.
//! * [`PackEngine`] — the runnable surface: the shared `with_*` builders and
//!   [`run`](PackEngine::run). `run` is a **required** method with no
//!   provided body, so there is exactly one lifecycle body in the crate
//!   (`Pipeline`'s) and no way for an engine to wrap itself in a second one.
//!   The three presets implement it in one line each.
//!
//! [`EngineSetup`] lives here because both its producer (the lifecycle) and
//! its only consumer ([`StageFactory::stages`]) do: one fact, one home.

use molrs::op::types::F;
use molrs::spatial::SimBox;

use crate::PackError;
use crate::Stage;
use crate::Target;
use crate::entry::setup::CellDecl;
use crate::entry::{PackSettings, State};
use crate::handler::{Handler, LogLevel};

/// Everything the lifecycle resolved before handing control to the stages:
/// the run's shared settings, the targets (post-broadcast), the space, and
/// the context shape. Borrowed — valid only inside [`StageFactory::stages`].
pub struct EngineSetup<'a> {
    /// The run's shared settings — the one ruler.
    ///
    /// A factory builds its stages from *these*, never from its own
    /// [`StageFactory::settings`]: a pipeline adopts the settings of the
    /// engine it was given (and refuses a second set from any other stage),
    /// so a factory's own copy is not the run's authority.
    pub settings: &'a PackSettings,
    pub targets: &'a [Target],
    pub cell: Option<SimBox>,
    pub maxmove_per_type: &'a [usize],
    pub ntype: usize,
    pub ntype_with_fixed: usize,
    pub ntotmol_free: usize,
    pub ntotat: usize,
    pub ntotat_free: usize,
}

/// A contributor of stages to a run.
///
/// The unit a [`Pipeline`](super::Pipeline) composes: it validates the
/// targets it will be given, states the shared settings it carries, hands
/// over its handlers, and builds its stages once the setup is resolved.
///
/// Implement this to plug a new algorithm into a pipeline without giving it
/// its own lifecycle. Implement [`PackEngine`] on top when it should also be
/// runnable on its own.
pub trait StageFactory {
    /// Reject targets this factory cannot handle — by name, never by
    /// silently switching to another algorithm. Called before any context is
    /// built.
    fn validate_targets(&self, _targets: &[Target]) -> Result<(), PackError> {
        Ok(())
    }

    /// The shared settings this factory carries.
    ///
    /// Inside a pipeline these are read only to *refuse* a second ruler: a
    /// factory handed to [`Pipeline::with_stage`](super::Pipeline::with_stage)
    /// with any knob off its default is rejected by name
    /// ([`PackError::PresetSettingsInsidePipeline`]).
    fn settings(&self) -> &PackSettings;

    /// Surrender the handlers this factory carries. Default: none.
    ///
    /// A pipeline **adopts** them — it never drops them — and they then
    /// observe the whole run, not just the stage that carried them in. The
    /// counterpart rule is the one above: handlers are adopted, non-default
    /// shared settings are refused by name.
    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        Vec::new()
    }

    /// Build this factory's stages for the resolved setup, in run order.
    ///
    /// Takes `&mut self` because a factory may hand its stages resources it
    /// owns (a seed's placements, bound optimizers); this is the one-shot
    /// handover, and the stages themselves stay re-entrant.
    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError>;
}

/// A boxed factory is a factory: `with_stage` accepts one that is already
/// boxed, and a pipeline stores its own that way.
impl<T: StageFactory + ?Sized> StageFactory for Box<T> {
    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        (**self).validate_targets(targets)
    }
    fn settings(&self) -> &PackSettings {
        (**self).settings()
    }
    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        (**self).take_handlers()
    }
    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        (**self).stages(setup)
    }
}

/// A runnable packing engine: the shared builders plus one verb.
///
/// [`run`](Self::run) is required, not provided — the lifecycle body lives in
/// exactly one place, [`Pipeline`](super::Pipeline), and every other
/// implementor delegates to it (`Pipeline::single(self).run(targets, n)`).
/// Consuming `self` makes an engine one shot by construction, which is what
/// makes its handler set impossible to lose silently.
pub trait PackEngine: StageFactory + Sized {
    /// Mutate the shared settings (used by the provided `with_*` builders).
    fn settings_mut(&mut self) -> &mut PackSettings;
    /// The engine's handler set.
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>>;

    // ── Shared builders (bound once, forwarded to `PackSettings`) ────────

    fn with_tolerance(mut self, tolerance: F) -> Self {
        self.settings_mut().tolerance = Some(tolerance);
        self
    }
    fn with_precision(mut self, precision: F) -> Self {
        self.settings_mut().precision = Some(precision);
        self
    }
    fn with_seed(mut self, seed: u64) -> Self {
        self.settings_mut().seed = Some(seed);
        self
    }
    fn with_periodic_box(mut self, min: [F; 3], max: [F; 3], pbc: [bool; 3]) -> Self {
        self.settings_mut().periodic_box = Some((min, max, pbc));
        self
    }
    fn with_density(mut self, rho: F) -> Self {
        self.settings_mut().density = Some(rho);
        self
    }
    fn with_parallel_eval(mut self, on: bool) -> Self {
        self.settings_mut().parallel_eval = on;
        self
    }
    fn with_log_level(mut self, level: LogLevel) -> Self {
        self.settings_mut().log.level = level;
        self
    }
    fn with_log_frequency(mut self, n: usize) -> Self {
        self.settings_mut().log.frequency = n.max(1);
        self
    }
    fn with_handler(mut self, handler: Box<dyn Handler>) -> Self {
        self.handlers_mut().push(handler);
        self
    }
    /// Declare the packing cell by lengths and angles (script `cell`).
    fn with_cell(mut self, lengths: [F; 3], angles_deg: [F; 3], pbc: [bool; 3]) -> Self {
        self.settings_mut().cell = Some(CellDecl::LengthsAngles {
            lengths,
            angles_deg,
            pbc,
        });
        self
    }
    /// Declare the packing cell by its lattice matrix.
    fn with_cell_matrix(mut self, h: [[F; 3]; 3], origin: [F; 3], pbc: [bool; 3]) -> Self {
        self.settings_mut().cell = Some(CellDecl::Matrix { h, origin, pbc });
        self
    }
    /// Packmol's optional second, shorter-range penalty. Declare the
    /// tolerance first: the distance must stay below it.
    fn with_short_tolerance(mut self, distance: F, scale: F) -> Self {
        assert!(
            distance > 0.0 && !distance.is_nan(),
            "short tolerance distance must be positive, got {distance}"
        );
        assert!(
            scale > 0.0 && !scale.is_nan(),
            "short tolerance scale must be positive, got {scale}"
        );
        let tolerance = self.settings().tolerance();
        assert!(
            distance < tolerance,
            "short tolerance distance {distance} must be smaller than the tolerance {tolerance}"
        );
        // Stored halved: the context wants the per-atom short radius.
        self.settings_mut().short_tolerance = Some((distance / 2.0, scale));
        self
    }
    /// Broadcast a restraint to every target at run time.
    fn with_global_restraint(mut self, r: impl crate::AtomRestraint + 'static) -> Self {
        self.settings_mut()
            .global_restraints
            .push(std::sync::Arc::new(r));
        self
    }

    /// Run the packing. Consumes the engine: one engine, one run.
    fn run(self, targets: &[Target], max_loops: usize) -> Result<State, PackError>;
}
