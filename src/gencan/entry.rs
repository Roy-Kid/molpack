//! `GenCanPack` — the rigid-body GENCAN packing entry.

use molrs::types::F;

use crate::entry::result::Placements;
use crate::entry::setup::CellDecl;
use crate::entry::{PackSettings, State};
use crate::error::PackError;
use crate::gencan::solver::{GencanSettings, GencanStage};
use crate::handler::{EarlyStopHandler, Handler};
use crate::optimizer::OptimizerBinding;
use crate::pipeline::{EngineSetup, PackEngine, Pipeline, StageFactory};
use crate::stage::Stage;
use crate::target::Target;

/// Rigid-body packing via the GENCAN bound-constrained optimizer
/// (Birgin & Martínez) — the Packmol algorithm as its own entry.
///
/// Shared knobs (`with_seed`, `with_tolerance`, boxes, handlers, …) come
/// from [`PackEngine`]; everything on this type is GENCAN-only and means
/// nothing to a growth entry.
pub struct GenCanPack {
    settings: PackSettings,
    handlers: Vec<Box<dyn Handler>>,
    early_stop: Option<EarlyStopHandler>,
    inner_iterations: usize,
    init_passes: Option<usize>,
    init_box_half_size: F,
    perturb_fraction: F,
    random_perturb: bool,
    perturb: bool,
    avoid_overlap: bool,
    seed_placements: Option<Placements>,
    optimizers: Vec<OptimizerBinding>,
}

impl Default for GenCanPack {
    fn default() -> Self {
        Self::new()
    }
}

impl GenCanPack {
    /// Packmol's default `nloop` for `ntype` structure types: `200 * ntype`
    /// GENCAN loops per phase (getinp.f90:537-539). What `max_loops` to pass
    /// when the caller has no reason to pick another.
    pub const fn default_max_loops(ntype: usize) -> usize {
        200 * ntype
    }

    pub fn new() -> Self {
        // One home for the GENCAN defaults: `GencanSettings::default()`.
        let d = GencanSettings::default();
        Self {
            settings: PackSettings::default(),
            handlers: Vec::new(),
            early_stop: Some(EarlyStopHandler::default()),
            inner_iterations: d.inner_iterations,
            init_passes: d.init_passes,
            init_box_half_size: d.init_box_half_size,
            perturb_fraction: d.perturb_fraction,
            random_perturb: d.random_perturb,
            perturb: d.perturb,
            avoid_overlap: d.avoid_overlap,
            seed_placements: None,
            optimizers: Vec::new(),
        }
    }

    /// Continue from a previous run's placement solution: the free copies
    /// start EXACTLY where `result` left them, bit for bit.
    ///
    /// Two things are switched off for such a run. The random-placement pass
    /// (`initial()`) is skipped, since the placements already exist; and so is
    /// movebad, the stall heuristic that picks up the worst-placed molecules
    /// and drops them elsewhere in the box. What is left is continuous
    /// rigid-body descent, which separates the remaining contacts while
    /// keeping the structure it was handed — the "slow push-off" of Auhl et al.
    /// 2003: an overlapped melt is opened up by continuous minimization, never
    /// by re-randomizing a configuration whose large-scale chain statistics
    /// were expensive to build.
    ///
    /// The cell travels with the seed; do not declare a box, density, or cell
    /// on a seeded engine. The run's free targets must match the seed's shape
    /// ([`PackError::SeedMismatch`] otherwise); fixed targets may be appended
    /// after the free ones.
    pub fn with_restart(mut self, result: &State) -> Self {
        // The seed's cell flows through the shared settings — one source of
        // truth, and the existing mutual-exclusion errors fire if the caller
        // declares a second box.
        self.settings.cell = Some(CellDecl::Resolved(result.placements.cell.clone()));
        self.seed_placements = Some(result.placements.clone());
        self
    }

    /// Bind an in-loop optimizer to a selection of copies.
    pub fn with_optimizer(
        mut self,
        select: crate::OptimizeSelect,
        optimizer: impl molrs::optimize::Optimizer + 'static,
    ) -> Self {
        self.optimizers.push(crate::optimizer::OptimizerBinding {
            select,
            optimizer: Box::new(optimizer),
        });
        self
    }

    /// GENCAN inner iterations (`maxit`).
    pub fn with_inner_iterations(mut self, n: usize) -> Self {
        self.inner_iterations = n.max(1);
        self
    }

    /// Initialization outer loops (`nloop0`); default `20 * ntype`.
    pub fn with_init_passes(mut self, n: usize) -> Self {
        self.init_passes = Some(n);
        self
    }

    /// Maximum system half-size in the initial restmol stage (`sidemax`).
    pub fn with_init_box_half_size(mut self, half_size: F) -> Self {
        self.init_box_half_size = half_size;
        self
    }

    /// Stall-perturbation heuristic: fraction perturbed / random selection /
    /// master switch (Packmol `movefrac` / `movebadrandom`).
    pub fn with_perturb(mut self, fraction: F, random: bool, enabled: bool) -> Self {
        self.perturb_fraction = fraction;
        self.random_perturb = random;
        self.perturb = enabled;
        self
    }

    /// End a phase whose best objective (Packmol's `bestf`) has stopped
    /// improving at the user's radii, instead of running it to `max_loops`.
    ///
    /// On by default with [`EarlyStopHandler::default`] — 10 % over 10 loops
    /// at `radscale == 1`, the criterion spelled out on that type. A pack
    /// that stalls there is almost always too dense, and the useful answer is
    /// `converged == false` in minutes, not an hour of GENCAN. Pass a
    /// configured handler to tune it, or `None` to run every phase to
    /// `max_loops` as Packmol does.
    pub fn with_early_stop(mut self, early_stop: impl Into<Option<EarlyStopHandler>>) -> Self {
        self.early_stop = early_stop.into();
        self
    }

    /// Reject initial placements that overlap a fixed molecule.
    pub fn with_avoid_overlap(mut self, on: bool) -> Self {
        self.avoid_overlap = on;
        self
    }
}

impl StageFactory for GenCanPack {
    fn settings(&self) -> &PackSettings {
        &self.settings
    }

    fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> {
        let Some(seed) = &self.seed_placements else {
            return Ok(());
        };
        // The seed covers the free copies; their shape must match per copy.
        let expected: Vec<usize> = targets
            .iter()
            .filter(|t| t.fixed_at.is_none())
            .flat_map(|t| std::iter::repeat_n(t.natoms(), t.count))
            .collect();
        if expected != seed.copy_atoms {
            return Err(PackError::SeedMismatch {
                expected: expected.iter().sum(),
                got: seed.copy_atoms.iter().sum(),
            });
        }
        Ok(())
    }

    fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> {
        let mut handlers = std::mem::take(self.handlers_mut());
        if let Some(early_stop) = self.early_stop.take() {
            handlers.push(Box::new(early_stop));
        }
        handlers
    }

    fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError> {
        // Shared knobs come from the run (`setup.settings`), algorithm knobs
        // from this entry — one ruler, one owner each.
        let s = setup.settings;
        let gencan = GencanSettings {
            inner_iterations: self.inner_iterations,
            init_passes: self.init_passes,
            init_box_half_size: self.init_box_half_size,
            perturb_fraction: self.perturb_fraction,
            random_perturb: self.random_perturb,
            perturb: self.perturb,
            avoid_overlap: self.avoid_overlap,
            discale: s.discale(),
            seed: s.seed(),
        };
        let stage = GencanStage::new(
            gencan,
            setup.maxmove_per_type.to_vec(),
            setup.cell.clone(),
            setup.ntype,
            setup.ntype_with_fixed,
        );
        // Both handovers are the entry's one-shot move, not a per-run one:
        // `run(self)` consumes the entry, so there is no second `stages()`
        // call to run bare, and the stage keeps what it is given for every
        // run it is asked to do.
        let stage = match self.seed_placements.take() {
            Some(seed) => stage.with_seed_placements(seed),
            None => stage,
        };
        let stage = stage.with_optimizers(std::mem::take(&mut self.optimizers));
        Ok(vec![Box::new(stage)])
    }
}

impl PackEngine for GenCanPack {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::target::Target;

    /// The entry lifecycle is deterministic: same targets, same seed,
    /// bit-identical positions and verdict. (Bitwise parity against the
    /// deleted `Molpack` path was proven before its removal —
    /// engine-entry-split migration record.)
    #[test]
    fn gencan_entry_is_deterministic() {
        let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let run = || {
            GenCanPack::new()
                .with_seed(11)
                .with_tolerance(2.0)
                .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
                .run(&[Target::from_coords(&coords, &[1.0, 1.0], 6)], 50)
                .expect("entry pack runs")
        };
        let (a, b) = (run(), run());
        assert!(a.converged);
        assert_eq!(a.fdist.to_bits(), b.fdist.to_bits());
        for (pa, pb) in a.positions().iter().zip(b.positions().iter()) {
            for k in 0..3 {
                assert_eq!(pa[k].to_bits(), pb[k].to_bits());
            }
        }
    }
}
