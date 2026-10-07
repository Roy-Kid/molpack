//! The GENCAN outer machinery: per-phase scaffold and per-iteration step.
//!
//! Free functions pulled out of the packer main loop (phases A.4.1-A.4.3)
//! and moved beside the optimizer they drive (engine-entry-split): the
//! stage owns the phase loop, these own one phase and one iteration.
//! Step reports use [`super::STAGE_NAME`], the same string the stage reports.

use molrs::op::F;
use rand::rngs::SmallRng;

use crate::Objective;
use crate::context::PackContext;
use crate::eval::EvalMode;
// The unscaled verdict is a shared primitive owned by the context layer, not
// by this stage: growth evaluates the same way, and the pipeline layer must
// not import `gencan/`.
use super::small_floor;
use crate::callback::{Callback, PhaseProgress, PhaseReport, StageProgress, StepReport};
use crate::context::pack_state::evaluate_unscaled;
use crate::optimizer::{ResolvedBinding, run_optimizer_bindings};
use crate::pack::gencan::{GencanParams, GencanWorkspace, pgencan};
use crate::pack::initial::SwapState;
use crate::pack::movebad::{MoveBadConfig, movebad};

/// Outcome of one main-loop iteration inside a packing phase.
///
/// The per-iteration body runs movebad → in-loop optimizers → pgencan → radii
/// schedule. `Continue`
/// means "run the next iteration"; `Converged` means the convergence predicate
/// fired inside this iteration; `EarlyStop` means a `Callback::should_stop()`
/// returned true.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IterOutcome {
    Continue,
    Converged,
    EarlyStop,
}

/// Run one iteration of a packing phase's main loop.
///
/// Matches Packmol's per-iteration sequence in `app/packmol.f90` lines 815-948:
///
/// 1. `movebad` when `radscale == 1.0` and previous `fimp <= 10%` (unless
///    disabled).
/// 2. Per-target in-loop optimizer block (`run_optimizer_bindings`).
/// 3. `pgencan` on the working coordinate vector.
/// 4. Unscaled-radii statistics (`fdist` / `frest` / `fimp`).
/// 5. Callback `on_step` notification; early stop if any callback opts in.
/// 6. Convergence check (`fdist < precision && frest < precision`).
/// 7. Radii reduction schedule (only when `radscale > 1.0`).
///
/// The function takes each piece of mutable outer-loop state by `&mut` so the
/// caller (the outer phase for-loop in `pack()`) retains ownership across
/// iterations.
#[allow(clippy::too_many_arguments)]
pub fn run_iteration(
    loop_idx: usize,
    max_loops: usize,
    is_all: bool,
    phase: usize,
    phase_progress: PhaseProgress,
    precision: F,
    disable_movebad: bool,
    movebad_cfg: &MoveBadConfig,
    gencan_params: &GencanParams,
    sys: &mut PackContext,
    xwork: &mut [F],
    swap: &mut SwapState,
    flast: &mut F,
    fimp_prev: &mut F,
    radscale: &mut F,
    optimizer_bindings: &mut [ResolvedBinding<'_>],
    callbacks: &mut [Box<dyn Callback>],
    gencan_workspace: &mut GencanWorkspace,
    rng: &mut SmallRng,
) -> IterOutcome {
    // movebad: Packmol triggers when radscale==1.0 AND fimp<=10.0
    // (packmol.f90 line 815). fimp here is from the PREVIOUS iteration.
    // After movebad, reset flast to the post-movebad f (Packmol line 821).
    if !disable_movebad && *radscale == 1.0 && *fimp_prev <= 10.0 {
        movebad(xwork, sys, precision, movebad_cfg, rng, gencan_workspace);
        // Reset flast to the post-movebad f value so fimp is measured
        // relative to movebad's starting point.
        *flast = evaluate_unscaled(sys, xwork).0;
    }

    // In-loop optimizers: all-type phase only (full x ⇒ clean COM/Euler).
    if is_all {
        run_optimizer_bindings(sys, xwork, optimizer_bindings);
    }

    // GENCAN on working x (compact for per-type, full for all-type)
    sys.reset_eval_counters();
    let res = pgencan(xwork, sys, gencan_params, precision, gencan_workspace);

    // Save compact results back to swap (for restore later)
    if !is_all {
        swap.save_type(phase, xwork, sys);
    }

    // Compute statistics with unscaled radii
    // (Packmol lines 833-841: radiuswork + computef + restore)
    let (fx_unscaled, fdist, frest) = evaluate_unscaled(sys, xwork);

    // fimp: percentage improvement in unscaled f from last iteration
    // Packmol line 846: if(flast>0) fimp = -100*(fx-flast)/flast
    let mut fimp = if *flast > 0.0 {
        -100.0 * (fx_unscaled - *flast) / *flast
    } else if fx_unscaled < small_floor() {
        100.0 // already converged
    } else {
        F::INFINITY
    };
    // Packmol lines 848-849: clamp to [-99.99, 99.99]
    fimp = fimp.clamp(-99.99, 99.99);
    *flast = fx_unscaled;
    *fimp_prev = fimp;

    if !callbacks.is_empty() {
        let step = StepReport {
            // One stage per run until the pipeline lands; the name is the
            // stage's own, taken from the stage type so the two cannot drift.
            stage: StageProgress {
                index: 0,
                total: 1,
                name: super::STAGE_NAME,
            },
            loop_idx,
            max_loops,
            phase: phase_progress,
            fdist,
            frest,
            f: fx_unscaled,
            improvement_pct: fimp,
            radscale: *radscale,
            precision,
        };
        for h in callbacks.iter_mut() {
            h.on_step(&step, sys);
        }

        if callbacks.iter().any(|h| h.should_stop()) {
            log::debug!("  Early stop requested at loop {loop_idx}");
            return IterOutcome::EarlyStop;
        }
    }

    log::debug!(
        "    loop={loop_idx} f={:.4e} fdist={:.4e} frest={:.4e} radscale={:.4} fimp={:.2}% ncf={} ncg={} inform={}",
        res.f,
        fdist,
        frest,
        *radscale,
        fimp,
        sys.ncf(),
        sys.ncg(),
        res.inform
    );

    // Check convergence
    if fdist < precision && frest < precision {
        log::debug!("  Converged at phase {phase} loop {loop_idx}");
        return IterOutcome::Converged;
    }

    // Radii reduction schedule (Packmol lines 940-948):
    //   if (fdist<precision && fimp<10%) || fimp<2%: reduce radscale
    if *radscale > 1.0 && (fimp < 2.0 || (fdist < precision && fimp < 10.0)) {
        *radscale = (0.9 * *radscale).max(1.0);
        for i in 0..sys.ntotat {
            let new_r = sys.radius_ini[i].max(0.9 * sys.radius[i]);
            sys.set_radius(i, new_r);
        }
    }

    IterOutcome::Continue
}

/// Outcome of one outer-loop phase.
///
/// The per-phase scaffold covers callback phase-start notification, comptype
/// reconfiguration, radii reset, swap setup, pre-loop precision
/// short-circuit, inner GENCAN loop, and swap restore / xwork-back copy.
/// `Continue` means the outer phase loop should
/// proceed; `Converged` means the all-type phase converged and the outer loop
/// should break.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseOutcome {
    Continue,
    Converged,
}

/// Run one phase of the main packing loop (per-type or all-type).
///
/// Matches the outer `for phase in 0..=ntype` body of Packmol `app/packmol.f90`
/// lines 740-990 (the swaptype / comptype dance bracketing the GENCAN inner
/// loop). For a per-type phase (`phase < ntype`), `xwork` is a compact
/// `nmols[phase] * 6`-element slice produced by `SwapState::set_type`; for the
/// all-type phase (`phase == ntype`), `xwork` is a full `6 * ntotmol_free`
/// clone of `x`.
///
/// The function takes the outer-loop state (`sys`, `x`, `swap`,
/// optimizer bindings, `callbacks`, `gencan_workspace`, `rng`) by `&mut` so that
/// state persists across phases, exactly as the inlined body did.
///
/// Returns `PhaseOutcome::Converged` **only** when the all-type phase
/// converges (either on its entry precision check or inside the inner loop);
/// every per-type phase returns `Continue` regardless of whether that type
/// converged on its own (Packmol lets the all-type phase decide).
#[allow(clippy::too_many_arguments)]
pub fn run_phase(
    phase: usize,
    ntype: usize,
    ntype_with_fixed: usize,
    total_phases: usize,
    max_loops: usize,
    discale: F,
    precision: F,
    disable_movebad: bool,
    movebad_cfg: &MoveBadConfig,
    gencan_params: &GencanParams,
    sys: &mut PackContext,
    x: &mut [F],
    swap: &mut SwapState,
    optimizer_bindings: &mut [ResolvedBinding<'_>],
    callbacks: &mut [Box<dyn Callback>],
    gencan_workspace: &mut GencanWorkspace,
    rng: &mut SmallRng,
) -> PhaseOutcome {
    let is_all = phase == ntype;

    let phase_progress = PhaseProgress {
        phase,
        total_phases,
        molecule_type: if is_all { None } else { Some(phase) },
    };

    // Reset callback state between phases (e.g. EarlyStopCallback stall counter)
    for h in callbacks.iter_mut() {
        h.on_phase_start(&phase_progress);
    }

    // Set comptype for this phase
    for itype in 0..ntype_with_fixed {
        sys.comptype[itype] = if is_all {
            true
        } else {
            itype >= ntype || itype == phase
        };
    }

    log::debug!(
        "  Packing phase {phase} ({})",
        if is_all {
            "all".to_string()
        } else {
            format!("type {phase}")
        }
    );

    // Compact x to this type (action=1) or restore full x (all-type phase)
    // Packmol resets radscale = discale at the START of each phase.
    let mut radscale = discale;
    for icart in 0..sys.ntotat {
        sys.set_radius(icart, discale * sys.radius_ini[icart]);
    }

    // Get working x vector (compact for per-type, full for all-type)
    let mut xwork: Vec<F> = if !is_all {
        // Compact: n = nmols[phase] * 6
        // Re-save current x (action=0) then compact (action=1)
        *swap = SwapState::init(x, sys);
        swap.set_type(phase, sys)
    } else {
        // All-type: restore full x (action=3), use it directly
        swap.restore(x, sys);
        x.to_vec()
    };

    // Packmol checks whether the current approximation is already a solution
    // before entering the GENCAN loop for this phase (packmol.f90 lines 775-782).
    sys.evaluate(&xwork, EvalMode::FOnly, None);
    if sys.fdist < precision && sys.frest < precision {
        let report = PhaseReport {
            iterations: 0,
            fdist: sys.fdist,
            frest: sys.frest,
            converged: true,
        };
        for h in callbacks.iter_mut() {
            h.on_phase_end(&phase_progress, &report);
        }
        if !is_all {
            swap.save_type(phase, &xwork, sys);
            swap.restore(x, sys);
            return PhaseOutcome::Continue;
        } else {
            x.copy_from_slice(&xwork);
            return PhaseOutcome::Converged;
        }
    }

    // Initialize flast = unscaled f before gencanloop
    // (Packmol lines 796-803: compute bestf/flast with unscaled radii)
    let mut flast = evaluate_unscaled(sys, &xwork).0;

    // fimp from previous iteration — used for movebad gate (Packmol packmol.f90 line 798).
    // Initialized to 1e99 so movebad is NOT called on the first iteration.
    let mut fimp_prev = F::INFINITY;
    let mut converged_inner = false;
    let mut iterations = 0usize;

    for loop_idx in 0..max_loops {
        let outcome = run_iteration(
            loop_idx,
            max_loops,
            is_all,
            phase,
            phase_progress,
            precision,
            disable_movebad,
            movebad_cfg,
            gencan_params,
            sys,
            &mut xwork,
            swap,
            &mut flast,
            &mut fimp_prev,
            &mut radscale,
            optimizer_bindings,
            callbacks,
            gencan_workspace,
            rng,
        );
        iterations += 1;
        match outcome {
            IterOutcome::Continue => {}
            IterOutcome::Converged => {
                converged_inner = true;
                break;
            }
            IterOutcome::EarlyStop => break,
        }
    }

    let report = PhaseReport {
        iterations,
        fdist: sys.fdist,
        frest: sys.frest,
        converged: converged_inner,
    };
    for h in callbacks.iter_mut() {
        h.on_phase_end(&phase_progress, &report);
    }

    // After per-type phase: save results + restore full x
    // After all-type phase: copy xwork back to x
    if !is_all {
        // save_type was called inside the loop; restore full x now.
        // Per-type convergence does NOT exit the outer phase loop.
        swap.restore(x, sys);
        PhaseOutcome::Continue
    } else {
        x.copy_from_slice(&xwork);
        if converged_inner {
            PhaseOutcome::Converged
        } else {
            PhaseOutcome::Continue
        }
    }
}
