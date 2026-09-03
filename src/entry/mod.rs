//! Per-algorithm engine entries and their shared lifecycle (engine-entry-split).
//!
//! One entry type per packing algorithm (`GenCanPack`, `CbmcGrow`, …), all
//! implementing [`PackEngine`]: the trait owns the lifecycle — validation,
//! space resolution, context construction, handler bracketing, assembly —
//! and each entry contributes exactly one thing, its [`Solver`].
//!
//! `run` consumes the entry by value: an engine is one shot by construction,
//! which is what makes the handler set impossible to lose silently (the old
//! `pack(&mut self)` drained its handlers on first use and ran headless on
//! the second).

pub(crate) mod result;
pub(crate) mod setup;

pub use result::PackResult;
pub(crate) use result::positions_in_target_order;

use molrs::spatial::simbox::SimBox;
use molrs::types::F;
use ndarray::Array1;

use crate::context::build::{ContextKnobs, build_context};
use crate::error::PackError;
use crate::handler::{Handler, LammpsLogHandler, LogLevel};
use crate::initial::init_xcart_from_x;

use crate::solver::{Budget, PlacementsMut, Solver};
use crate::target::Target;
use setup::{CellDecl, PeriodicSpec, broadcast_global_restraints, resolve_pack_space};

/// Built-in screen logging: detail level + print cadence.
#[derive(Debug, Clone, Copy)]
pub struct LogSpec {
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

/// The knobs every entry shares — consumed by the lifecycle and the shared
/// infrastructure (space resolution, context construction), never by one
/// algorithm alone. Algorithm-specific knobs live on their entry.
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
    pub(crate) global_restraints: Vec<std::sync::Arc<dyn crate::restraint::AtomRestraint>>,
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
}

/// Everything the lifecycle resolved before handing control to the entry's
/// [`Solver`]: the targets (post-broadcast), the space, and the context
/// shape. Borrowed — valid only inside [`PackEngine::solver`] /
/// [`PackEngine::prepare`].
pub struct EngineSetup<'a> {
    pub targets: &'a [Target],
    pub cell: Option<SimBox>,
    pub maxmove_per_type: &'a [usize],
    pub ntype: usize,
    pub ntype_with_fixed: usize,
    pub ntotmol_free: usize,
    pub ntotat: usize,
    pub ntotat_free: usize,
}

/// One packing algorithm behind one lifecycle.
///
/// Entries supply their [`Solver`] (and optionally validation and context
/// preparation); the provided [`run`](Self::run) owns everything shared:
/// space resolution, restraint broadcast, context construction, handler
/// bracketing (`on_start` / `on_finish`, log-handler injection), and frame
/// assembly. `fdist` / `frest` come from the shared objective — the seam's
/// one-ruler guarantee.
pub trait PackEngine: Sized {
    /// Read the shared settings.
    fn settings(&self) -> &PackSettings;
    /// Mutate the shared settings (used by the provided `with_*` builders).
    fn settings_mut(&mut self) -> &mut PackSettings;
    /// The entry's handler set.
    fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>>;
    /// Entry-specific target validation — named rejections only, never a
    /// silent fallback to another algorithm.
    fn validate(&self, _targets: &[Target]) -> Result<(), PackError> {
        Ok(())
    }
    /// Entry-specific context preparation (e.g. growth installs the box and
    /// cell grid before its solver runs; a seeded GENCAN run injects the
    /// seed placements into `sys.coor` / `x`; GENCAN's `initial()` does its
    /// own otherwise).
    fn prepare(
        &self,
        _sys: &mut crate::context::PackContext,
        _x: &mut [F],
        _setup: &EngineSetup<'_>,
    ) -> Result<(), PackError> {
        Ok(())
    }
    /// The algorithm: build this entry's [`Solver`] for the resolved setup.
    fn solver(&mut self, setup: &EngineSetup<'_>) -> Result<Box<dyn Solver>, PackError>;

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
    fn with_global_restraint(mut self, r: impl crate::restraint::AtomRestraint + 'static) -> Self {
        self.settings_mut()
            .global_restraints
            .push(std::sync::Arc::new(r));
        self
    }

    /// Run the packing. Consumes the entry: one engine, one run.
    fn run(mut self, targets: &[Target], max_loops: usize) -> Result<PackResult, PackError> {
        if targets.is_empty() {
            return Err(PackError::NoTargets);
        }
        for (i, t) in targets.iter().enumerate() {
            if t.natoms() == 0 {
                return Err(PackError::EmptyMolecule(i));
            }
        }
        self.validate(targets)?;

        // Stage ①: space + restraint broadcast (shared).
        let settings = self.settings();
        let broadcast = broadcast_global_restraints(targets, &settings.global_restraints);
        let targets: &[Target] = &broadcast;
        let space = resolve_pack_space(
            settings.density,
            settings.periodic_box,
            settings.cell,
            targets,
        )?;

        // Stage ②: context.
        let knobs = ContextKnobs {
            tolerance: settings.tolerance(),
            short_tolerance: settings.short_tolerance,
            parallel_eval: settings.parallel_eval,
        };
        let built = build_context(&knobs, targets)?;
        let mut sys = built.sys;
        let mut x = vec![0.0 as F; 6 * built.ntotmol_free];

        let setup = EngineSetup {
            targets,
            cell: space.cell.clone(),
            maxmove_per_type: &built.maxmove_per_type,
            ntype: built.ntype,
            ntype_with_fixed: built.ntype_with_fixed,
            ntotmol_free: built.ntotmol_free,
            ntotat: built.ntotat,
            ntotat_free: built.ntotat_free,
        };

        // Stages ③④: the algorithm, bracketed by the handlers.
        let mut solver = self.solver(&setup)?;
        self.prepare(&mut sys, &mut x, &setup)?;

        let log = self.settings().log;
        let (tolerance, precision, seed) = {
            let s = self.settings();
            (s.tolerance(), s.precision(), s.seed())
        };
        let mut handlers = std::mem::take(self.handlers_mut());
        if log.level.is_enabled() {
            handlers.push(Box::new(LammpsLogHandler::new(
                log.level,
                log.frequency,
                tolerance,
                precision,
                seed,
                max_loops,
                built.ntype_with_fixed,
                space.pbc,
            )));
        }
        for h in handlers.iter_mut() {
            h.on_start(built.ntotat, built.ntotmol_free);
        }
        sys.ntotmol = built.ntotmol_free;

        let budget = Budget::new(max_loops, precision);
        let outcome = solver.solve(
            &mut sys,
            targets,
            PlacementsMut::new(&mut x, built.ntotmol_free),
            &budget,
            &mut handlers,
        );
        // Stage ⑤: xcart rebuild, finish bracket, assembly.
        for itype in 0..built.ntype_with_fixed {
            sys.comptype[itype] = true;
        }
        sys.ntotmol = built.ntotmol_free;
        init_xcart_from_x(&x, &mut sys);
        for h in handlers.iter_mut() {
            h.on_finish(&sys);
        }

        // The placement solution, captured verbatim for cross-entry seeding
        // (placement-seeding spec): the frame below is derived VIEW data —
        // re-deriving (coor, x) from it would recompute COMs and break
        // bitwise continuity.
        let placements = result::Placements {
            x: x.clone(),
            coor: sys.coor[..built.ntotat_free].to_vec(),
            copy_atoms: targets
                .iter()
                .filter(|t| t.fixed_at.is_none())
                .flat_map(|t| std::iter::repeat_n(t.natoms(), t.count))
                .collect(),
            cell: sys.simbox.clone(),
        };

        let xcart = std::mem::take(&mut sys.xcart);
        let positions = positions_in_target_order(targets, &xcart, built.ntotat_free);
        let mut frame = crate::assemble::assemble_frame(targets, &positions);
        if let Some((min, max, flags)) = space.pbc {
            let lengths = Array1::from_vec(vec![max[0] - min[0], max[1] - min[1], max[2] - min[2]]);
            let origin = Array1::from_vec(min.to_vec());
            if let Ok(simbox) = SimBox::ortho(lengths, origin, flags) {
                frame.simbox = Some(simbox);
            }
        }
        // A DECLARED cell (with_cell / a seeded run's inherited cell) is
        // user-stated geometry and belongs on the output frame; a box merely
        // inferred from restraints stays off it, as documented.
        if frame.simbox.is_none()
            && let Some(cell) = &space.cell
        {
            frame.simbox = Some(cell.clone());
        }

        Ok(PackResult {
            frame,
            placements,
            fdist: sys.fdist,
            frest: sys.frest,
            converged: outcome.converged,
            softened: outcome.softened,
        })
    }
}
