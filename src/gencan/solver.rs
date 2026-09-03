//! GENCAN on the [`Solver`] seam.
//!
//! Until engine-entry-split, the rigid-body path was a 12-argument inherent
//! method on the `Molpack` builder (`run_gencan_stages`) — the peer claim in
//! [`solver`](crate::solver) existed only in prose. [`GencanSolver`] makes it
//! true in code: the same lifecycle the growth solver implements, judged by
//! the same shared-objective ruler, selected by the same seam.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::context::{PackContext, RigidView};
use crate::gencan::phases::{PhaseOutcome, run_phase};
use crate::gencan::{GencanParams, GencanWorkspace};
use crate::handler::Handler;
use crate::initial::{SwapState, initial};
use crate::movebad::MoveBadConfig;
#[cfg(feature = "ff")]
use crate::optimizer::{OptimizerBinding, ResolvedBinding, resolve_bindings};
use crate::solver::{Budget, SolveOutcome, Solver};
use crate::target::Target;

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
    /// Continue from a pre-seeded state (skip `initial()` and movebad —
    /// the coordinates are already placed, e.g. by a growth stage). Not a
    /// `GenCanPack` builder knob; set by chaining code.
    pub push_off: bool,
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
            push_off: false,
        }
    }
}

/// The rigid-body GENCAN packing algorithm as a [`Solver`].
///
/// Construction captures everything the former `run_gencan_stages` read
/// beyond the seam signature: the GENCAN knobs, the per-type move quota,
/// the resolved cell, the phase-shape counts, and the push-off flag (skip
/// `initial()` and movebad — the entry seeded the state from a previous run,
/// see [`GenCanPack::seeded_from`](crate::GenCanPack::seeded_from)).
pub struct GencanSolver {
    settings: GencanSettings,
    maxmove_per_type: Vec<usize>,
    cell: Option<SimBox>,
    ntype: usize,
    ntype_with_fixed: usize,
    push_off: bool,
    #[cfg(feature = "ff")]
    optimizers: Vec<OptimizerBinding>,
    rng: SmallRng,
}

impl GencanSolver {
    pub fn new(
        settings: GencanSettings,
        maxmove_per_type: Vec<usize>,
        cell: Option<SimBox>,
        ntype: usize,
        ntype_with_fixed: usize,
    ) -> Self {
        let rng = SmallRng::seed_from_u64(settings.seed);
        let push_off = settings.push_off;
        Self {
            settings,
            maxmove_per_type,
            cell,
            ntype,
            ntype_with_fixed,
            push_off,
            #[cfg(feature = "ff")]
            optimizers: Vec::new(),
            rng,
        }
    }

    /// Bind in-loop optimizers (feature `ff`). Kept off the constructor so
    /// non-`gencan` callers never spell the `ff` cfg (acceptance gate:
    /// `src/grow/` stays free of `cfg(feature = "ff")`).
    #[cfg(feature = "ff")]
    pub fn with_optimizers(mut self, optimizers: Vec<OptimizerBinding>) -> Self {
        self.optimizers = optimizers;
        self
    }
}

impl Solver for GencanSolver {
    fn name(&self) -> &'static str {
        "gencan"
    }

    fn solve(
        &mut self,
        sys: &mut PackContext,
        targets: &[Target],
        x: &mut RigidView,
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> SolveOutcome {
        let init_passes = self.settings.init_passes.unwrap_or(20 * self.ntype);
        let movebad_cfg = MoveBadConfig {
            movefrac: self.settings.perturb_fraction,
            maxmove_per_type: &self.maxmove_per_type,
            movebadrandom: self.settings.random_perturb,
            gencan_maxit: self.settings.inner_iterations,
        };
        if !self.push_off {
            initial(
                x.as_mut_slice(),
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

        #[cfg(feature = "ff")]
        let mut optimizer_bindings: Vec<ResolvedBinding> = {
            let type_names: Vec<Option<String>> = targets
                .iter()
                .filter(|t| t.fixed_at.is_none())
                .map(|t| t.name.clone())
                .collect();
            resolve_bindings(std::mem::take(&mut self.optimizers), &type_names)
        };
        #[cfg(not(feature = "ff"))]
        let _ = targets;

        // max_loops controls the outer loop count, matching Packmol's `nloop`.
        let gencan_params = GencanParams {
            maxit: self.settings.inner_iterations,
            maxfc: self.settings.inner_iterations * 10,
            iprint: 0,
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
                self.push_off || !self.settings.perturb,
                &movebad_cfg,
                &gencan_params,
                sys,
                x.as_mut_slice(),
                &mut swap,
                #[cfg(feature = "ff")]
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

        SolveOutcome::new(converged, sys.fdist, sys.frest, 0)
    }
}

#[cfg(test)]
mod tests {
    use molrs::spatial::simbox::SimBox;
    use ndarray::Array1;

    use super::*;
    use crate::context::build::{ContextKnobs, build_context};

    /// RED-1 (engine-entry-split): the rigid-body path must run behind the
    /// `Solver` seam — same context plumbing as any other solver, verdict
    /// from the shared objective, no `Molpack` internals.
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
        let mut sys = built.sys;
        let mut x = RigidView::fresh(built.ntotmol_free);

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
        let mut solver: Box<dyn Solver> = Box::new(GencanSolver::new(
            settings,
            built.maxmove_per_type.clone(),
            Some(cell),
            built.ntype,
            built.ntype_with_fixed,
        ));

        let mut handlers: Vec<Box<dyn Handler>> = Vec::new();
        let outcome = solver.solve(
            &mut sys,
            &targets,
            &mut x,
            &Budget::new(50, 0.01),
            &mut handlers,
        );

        assert_eq!(solver.name(), "gencan");
        assert!(outcome.converged, "6 dimers in a 20 Å box must converge");
        assert_eq!(outcome.softened, 0, "GENCAN never softens");
        assert!(
            sys.fdist <= 0.01,
            "verdict comes from the shared objective: fdist = {}",
            sys.fdist
        );
    }
}
