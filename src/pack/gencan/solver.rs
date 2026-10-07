//! GENCAN on the [`Stage`] seam.
//!
//! [`GenCanStage`] is the rigid-body path behind the stage seam: the same
//! lifecycle the growth stages implement, judged by
//! the same shared-objective ruler, selected by the same seam.

use molrs::core::SimBox;
use molrs::op::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::Handler;
use crate::PackError;
use crate::Target;
use crate::context::{PackState, Placed, RigidView};
use crate::entry::result::Placements;
use crate::optimizer::{OptimizerBinding, ResolvedBinding, resolve_bindings};
use crate::pack::gencan::phases::{PhaseOutcome, run_phase};
use crate::pack::gencan::{GencanParams, GencanWorkspace};
use crate::pack::initial::{SwapState, initial};
use crate::pack::movebad::MoveBadConfig;
use crate::stage::{Budget, Guarantees, Requires, Stage, StageOutcome};

/// GENCAN-only knobs (engine-entry-split: these live on `GenCanPack`, never
/// on the shared settings — they mean nothing to a growth entry).
#[derive(Debug, Clone)]
pub struct GencanSettings {
    /// GENCAN inner iterations (`maxit`).
    pub inner_iterations: usize,
    /// Initialization outer loops (`nloop0`); `None` = Packmol's `20 * ntype`.
    pub init_passes: Option<usize>,
    /// Maximum system half-size in the initial restmol stage (`sidemax`).
    pub init_box_half_size: F,
    /// Fraction of molecules perturbed when packing stalls (`movefrac`).
    pub perturb_fraction: F,
    /// Randomize perturbation target selection (`movebadrandom`).
    pub random_perturb: bool,
    /// Master switch for the stall-perturbation heuristic.
    pub perturb: bool,
    /// Reject initial placements overlapping a fixed molecule (`avoid_overlap`).
    pub avoid_overlap: bool,
    /// Initial radius up-scaling (`discale`).
    pub discale: F,
    /// RNG seed; the solver owns its stream (bit-parity with the packer's
    /// former single `SmallRng`, which reached the GENCAN stage undrawn).
    pub seed: u64,
}

impl Default for GencanSettings {
    fn default() -> Self {
        Self {
            inner_iterations: 20,
            init_passes: None,
            init_box_half_size: 1000.0,
            perturb_fraction: 0.05,
            random_perturb: false,
            perturb: true,
            avoid_overlap: true,
            discale: 1.1,
            seed: 1_234_567,
        }
    }
}

/// The rigid-body GENCAN packing algorithm as a [`Stage`].
///
/// Construction captures everything the former `run_gencan_stages` read
/// beyond the seam signature: the GENCAN knobs, the per-type move quota,
/// the resolved cell, and the phase-shape counts. Whether the stage starts
/// from scratch or continues from placements it was handed is **not** stored
/// here — it is read off the state on entry (see [`run`](Stage::run)), which
/// is the only home that fact has.
///
/// The optimizer bindings and the placement seed stay on the stage for its
/// whole life and are read afresh on every [`run`](Stage::run) — the seam's
/// re-entrancy contract: a stage does not consume its own configuration.
pub struct GenCanStage {
    settings: GencanSettings,
    maxmove_per_type: Vec<usize>,
    cell: Option<SimBox>,
    ntype: usize,
    ntype_with_fixed: usize,
    seed_placements: Option<Placements>,
    optimizers: Vec<OptimizerBinding>,
    rng: SmallRng,
}

impl GenCanStage {
    /// The name this stage reports. Same constant as [`super::STAGE_NAME`],
    /// which the phase step report fills in, so the two cannot drift apart.
    pub(crate) const NAME: &'static str = super::STAGE_NAME;

    pub fn new(
        settings: GencanSettings,
        maxmove_per_type: Vec<usize>,
        cell: Option<SimBox>,
        ntype: usize,
        ntype_with_fixed: usize,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(settings.seed);
        Self {
            settings,
            maxmove_per_type,
            cell,
            ntype,
            ntype_with_fixed,
            seed_placements: None,
            optimizers: Vec::new(),
            rng,
        }
    }

    /// Continue from a previous run's placement solution.
    ///
    /// Kept off the constructor for the same reason as the optimizer
    /// bindings: it is an option of the seeded spelling only, and the
    /// stage is perfectly usable without it. Installed by
    /// [`run`](Stage::run), not here, so a second run re-injects the same
    /// bits rather than silently running degraded.
    pub(crate) fn with_seed_placements(mut self, placements: Placements) -> Self {
        self.seed_placements = Some(placements);
        self
    }

    /// Bind in-loop optimizers. Kept off the constructor like the placement
    /// seed: an option of the optimizer-bound spelling only.
    pub fn with_optimizers(mut self, optimizers: Vec<OptimizerBinding>) -> Self {
        self.optimizers = optimizers;
        self
    }
}

impl Stage for GenCanStage {
    fn name(&self) -> &'static str {
        Self::NAME
    }

    /// Nothing: GENCAN can place from scratch with `initial()`. When the
    /// state *does* arrive placed, it continues from those placements
    /// instead — see the preamble on [`run`](Stage::run).
    fn requires(&self) -> Requires {
        Requires::new(Placed::None)
    }

    /// Every free molecule placed.
    fn guarantees(&self) -> Guarantees {
        Guarantees::new(Placed::All)
    }

    /// # Preamble
    ///
    /// Three things happen before the packing loop, all of them decided in
    /// state vocabulary rather than by a flag the caller set:
    ///
    /// 1. **The box and the cell grid** are installed only when this run will
    ///    *continue* from placements that already exist — the state says
    ///    [`Placed::All`], or this stage carries a placement seed (which makes
    ///    it so in step 2). Starting from nothing, `initial()` owns the box
    ///    and the grid itself (synthesizing a fall-back box from `sidemax`
    ///    when nothing was declared), so installing one here would be a second
    ///    owner. The box is the entry's resolved cell when there is one, else
    ///    the one the context already carries — which is how a GENCAN stage
    ///    that follows another stage in a chain lands on the same box, and the
    ///    same `radmax`, as the hand-written `with_restart` spelling.
    /// 2. **The seed**, if this stage carries one, is injected verbatim
    ///    ([`RigidView::install_seed`](crate::RigidView::install_seed)) and the
    ///    state's marker advances to [`Placed::All`]. Order matters: the grid
    ///    first, the coordinates second. The seed is borrowed, never taken —
    ///    a second run re-injects the same bits (the seam's re-entrancy
    ///    contract).
    /// 3. **Push-off** is then simply `state.placed() == Placed::All`: skip
    ///    `initial()`, materialize `xcart` from the placements already there,
    ///    and keep movebad off so molecules are separated by descent rather
    ///    than teleported. No flag records it anywhere: the state is the one
    ///    home of that fact.
    fn run(
        &mut self,
        state: &mut PackState,
        targets: &[Target],
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> Result<StageOutcome, PackError> {
        // ① Box + cell grid, only for a run that continues from placements.
        if state.placed() == Placed::All || self.seed_placements.is_some() {
            let sys = state.ctx_mut();
            let simbox = self.cell.clone().unwrap_or_else(|| sys.simbox.clone());
            // One derivation of the grid's coverage scale, in
            // `initial::coverage_radmax`.
            crate::context::grid::install_resolved_cell(sys, &simbox, self.settings.discale);
        }
        // ② The seed's conformers and placements, verbatim.
        if let Some(seed) = &self.seed_placements {
            let (sys, x) = state.rigid_split_mut();
            *x = RigidView::install_seed(seed.rigid.as_slice(), &seed.coor, sys);
            state.set_placed(Placed::All);
        }
        // ③ Continue from what is placed, or place from nothing.
        let push_off = state.placed() == Placed::All;

        let (sys, x) = state.rigid_split_mut();
        let init_passes = self.settings.init_passes.unwrap_or(20 * self.ntype);
        let movebad_cfg = MoveBadConfig {
            movefrac: self.settings.perturb_fraction,
            maxmove_per_type: &self.maxmove_per_type,
            movebadrandom: self.settings.random_perturb,
            gencan_maxit: self.settings.inner_iterations,
        };
        if !push_off {
            initial(
                x,
                sys,
                budget.precision,
                self.settings.discale,
                self.settings.init_box_half_size,
                init_passes,
                self.cell.clone(),
                self.settings.avoid_overlap,
                &movebad_cfg,
                &mut self.rng,
            );
        } else {
            // Push-off: `x` and `sys.coor` were seeded from a previous run
            // and the grid is installed — materialize xcart from that state
            // and keep it; movebad stays disabled so molecules are pushed
            // apart by descent, never teleported.
            x.write_xcart(sys);
        }

        // Notify handlers: initialization complete, xcart is valid
        for h in handlers.iter_mut() {
            h.on_initialized(sys);
        }

        let mut optimizer_bindings: Vec<ResolvedBinding<'_>> = {
            let type_names: Vec<Option<String>> = targets
                .iter()
                .filter(|t| t.fixed_at.is_none())
                .map(|t| t.name.clone())
                .collect();
            // Borrowed, never taken: the bindings are this stage's own
            // configuration and the next run must find them intact.
            resolve_bindings(&mut self.optimizers, &type_names)
        };

        // max_loops controls the outer loop count, matching Packmol's `nloop`.
        let gencan_params = GencanParams {
            maxit: self.settings.inner_iterations,
            maxfc: self.settings.inner_iterations * 10,
            ..Default::default()
        };

        let mut converged = false;
        let mut gencan_workspace = GencanWorkspace::new();

        // ── Main optimization loop ─────────────────────────────────────────
        //
        // Matches Packmol's `app/packmol.f90` main loop exactly:
        //   For each type (itype 1..ntype): swaptype(action=1) → pack → restore
        //   Then all types (itype = ntype+1): pack with full x
        //
        // Per-type phases use a compact x (n = nmols[itype]*6) via SwapState,
        // reducing GENCAN problem size by up to 60x vs full n.

        // Save initial full x before phasing (Packmol swaptype action=0)
        let mut swap = SwapState::init(x.as_slice(), sys);

        let total_phases = self.ntype + 1;

        for phase in 0..=(self.ntype) {
            let outcome = run_phase(
                phase,
                self.ntype,
                self.ntype_with_fixed,
                total_phases,
                budget.max_loops,
                self.settings.discale,
                budget.precision,
                push_off || !self.settings.perturb,
                &movebad_cfg,
                &gencan_params,
                sys,
                x.as_mut_slice(),
                &mut swap,
                &mut optimizer_bindings,
                handlers,
                &mut gencan_workspace,
                &mut self.rng,
            );
            match outcome {
                PhaseOutcome::Continue => {}
                PhaseOutcome::Converged => {
                    converged = true;
                    break;
                }
            }
        }

        if !converged {
            log::warn!(
                "  Pack did not fully converge (fdist={:.4e}, frest={:.4e})",
                sys.fdist,
                sys.frest
            );
        }

        Ok(StageOutcome::new(converged, 0))
    }
}

#[cfg(test)]
mod tests {
    use molrs::core::SimBox;
    use ndarray::Array1;

    use super::*;
    use crate::context::build::{ContextKnobs, build_context};

    /// RED-1 (engine-entry-split): the rigid-body path must run behind the
    /// `Stage` seam — same context plumbing as any other stage, verdict
    /// from the shared objective, no entry internals.
    #[test]
    fn gencan_solves_a_small_pack_on_the_seam() {
        let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let targets = vec![Target::from_coords(&coords, &[1.0, 1.0], 6)];

        let built = build_context(
            &ContextKnobs {
                tolerance: 2.0,
                short_tolerance: None,
                parallel_eval: false,
            },
            &targets,
        )
        .expect("context builds");
        let mut state = PackState::new(built.sys, built.ntotmol_free);

        let cell = SimBox::ortho(
            Array1::from_vec(vec![20.0, 20.0, 20.0]),
            Array1::from_vec(vec![0.0, 0.0, 0.0]),
            [true; 3],
        )
        .expect("orthorhombic box");

        let settings = GencanSettings {
            seed: 7,
            ..Default::default()
        };
        let mut stage: Box<dyn Stage> = Box::new(GenCanStage::new(
            settings,
            built.maxmove_per_type.clone(),
            Some(cell),
            built.ntype,
            built.ntype_with_fixed,
        ));

        let mut handlers: Vec<Box<dyn Handler>> = Vec::new();
        let outcome = stage
            .run(&mut state, &targets, &Budget::new(50, 0.01), &mut handlers)
            .expect("gencan stage runs");

        assert_eq!(stage.name(), "gencan");
        assert!(outcome.converged, "6 dimers in a 20 Å box must converge");
        assert_eq!(outcome.degraded, 0, "GENCAN never softens");
        assert!(
            state.ctx().fdist <= 0.01,
            "verdict comes from the shared objective: fdist = {}",
            state.ctx().fdist
        );
    }

    /// ac-005 (stage-pipeline-04-stage): a stage may be run more than once on
    /// an evolving state, and the second run must have the same capabilities
    /// as the first — it may consume the scratch it builds per run, never its
    /// own configuration.
    ///
    /// **RED for the right reason.** Before this spec, `solve` resolved its
    /// bindings with `resolve_bindings(std::mem::take(&mut self.optimizers),
    /// ..)` (`src/gencan/solver.rs:176`), which *moves* the bindings off the
    /// stage. From the second run on, `self.optimizers` is empty, the
    /// optimizer block is silently a no-op, and the counter below stops
    /// advancing — a degraded result with no name (law § 10). The fix is to
    /// borrow the bindings rather than take them, after which the second run
    /// calls the optimizer exactly as the first did.
    ///
    /// **Why it lives in the crate.** `build_context` is `pub(crate)`, so an
    /// integration test in `tests/` cannot build a context and run the same
    /// stage twice on it; `optimizer::torsion_mc`'s tests only reaches the entry, which
    /// runs a stage once. Named with `optimizer` so the acceptance filter
    /// `cargo test -p molcrafts-molpack --lib -- optimizer` selects it.
    ///
    /// **Why the fixture is unsatisfiable.** The in-loop optimizer block runs
    /// only inside the all-type phase's iteration loop
    /// (`run_iteration` in `src/pack/gencan/phases.rs`), and a phase that is already a solution
    /// short-circuits past that loop. Twelve unit-radius dimers restrained
    /// into a 4 Å cube cannot be solved, so the loop is entered on both runs.
    /// The `after_first > 0` assertion guards exactly that: if the fixture
    /// ever stops reaching the optimizer, this test says so instead of
    /// certifying nothing.
    #[test]
    fn optimizer_bindings_survive_a_second_run() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicUsize, Ordering};

        use crate::OptimizeSelect;
        use crate::PackState;
        use crate::restraint::geometric::InsideBoxRestraint;

        /// An optimizer that only counts its calls and leaves the frame
        /// alone: the conformer is unchanged, so the non-harm gate is a
        /// no-op and the counter is the single thing this test observes.
        struct CountingOptimizer {
            calls: Arc<AtomicUsize>,
        }

        impl molrs::optimize::Optimizer for CountingOptimizer {
            fn minimize(
                &mut self,
                _frame: &mut molrs::core::Frame,
            ) -> Result<molrs::optimize::OptimizationReport, String> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                Ok(molrs::optimize::OptimizationReport {
                    converged: true,
                    n_steps: 0,
                    final_energy: 0.0,
                    final_fmax: 0.0,
                    final_grad_rms: 0.0,
                })
            }
        }

        let coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let targets = vec![
            Target::from_coords(&coords, &[1.0, 1.0], 12)
                .with_name("dimer")
                .with_restraint(InsideBoxRestraint::new([0.0; 3], [4.0; 3])),
        ];

        let built = build_context(
            &ContextKnobs {
                tolerance: 2.0,
                short_tolerance: None,
                parallel_eval: false,
            },
            &targets,
        )
        .expect("context builds");

        let cell = SimBox::ortho(
            Array1::from_vec(vec![20.0, 20.0, 20.0]),
            Array1::from_vec(vec![0.0, 0.0, 0.0]),
            [true; 3],
        )
        .expect("orthorhombic box");

        let calls = Arc::new(AtomicUsize::new(0));
        let settings = GencanSettings {
            seed: 7,
            ..Default::default()
        };
        let mut stage = GenCanStage::new(
            settings,
            built.maxmove_per_type.clone(),
            Some(cell),
            built.ntype,
            built.ntype_with_fixed,
        )
        .with_optimizers(vec![OptimizerBinding {
            select: OptimizeSelect::per_copy(["dimer"]),
            optimizer: Box::new(CountingOptimizer {
                calls: Arc::clone(&calls),
            }),
        }]);

        let mut state = PackState::new(built.sys, built.ntotmol_free);
        let budget = Budget::new(2, 0.01);
        let mut handlers: Vec<Box<dyn Handler>> = Vec::new();

        stage
            .run(&mut state, &targets, &budget, &mut handlers)
            .expect("gencan stage runs");
        let after_first = calls.load(Ordering::Relaxed);
        assert!(
            after_first > 0,
            "fixture guard: the first run never reached the in-loop optimizer, \
             so this test cannot say anything about the second one"
        );

        stage
            .run(&mut state, &targets, &budget, &mut handlers)
            .expect("gencan stage runs a second time");
        let after_second = calls.load(Ordering::Relaxed);

        assert!(
            after_second > after_first,
            "the optimizer counter did not advance during the second run \
             (after run 1: {after_first}, after run 2: {after_second}) — a \
             stage that consumes its own optimizer bindings on run 1 \
             degrades silently on every later run"
        );
    }
}
