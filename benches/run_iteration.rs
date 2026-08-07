//! Regression microbench for `packer::run_iteration` — one outer-loop step.
//!
//! `run_iteration` is the packer's per-iteration sequence (movebad → in-loop
//! optimizers → pgencan → unscaled statistics → handler notify → convergence
//! check → radii schedule). This bench drives one such step so a regression in the
//! iteration scaffold's boundary cost (indirection, inlining, the GENCAN
//! step) shows up.
//!
//! Setup: an empty-molecule `PackContext` (`ntotat = 4`, `ntotmol = 0`). With
//! `ntotmol = 0` the pgencan / evaluate / movebad branches run on empty
//! vectors, so this measures function-call boundary cost on a trivial body,
//! not full-workload throughput (that lives in `pack_end_to_end`). Fresh
//! snapshot per iteration via `iter_batched`, since a step mutates its state.
//!
//! Tiny and self-contained: no external test data, no features.

use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use molpack::gencan::{GencanParams, GencanWorkspace};
use molpack::handler::{Handler, PhaseInfo};
use molpack::initial::SwapState;
use molpack::movebad::MoveBadConfig;
use molpack::packer::{IterOutcome, run_iteration};
use molpack::{F, PackContext};

type Snapshot = (
    PackContext,
    Vec<F>,
    SwapState,
    GencanWorkspace,
    Vec<Box<dyn Handler>>,
    SmallRng,
);

fn build_snapshot() -> Snapshot {
    let ntotat = 4;
    let mut sys = PackContext::new(ntotat, 0, 0);
    sys.radius.fill(0.75);
    sys.radius_ini.fill(1.5);
    sys.work.radiuswork.resize(ntotat, 0.0);
    sys.sync_atom_props();
    let x: Vec<F> = Vec::new();
    let swap = SwapState::init(&x, &sys);
    let ws = GencanWorkspace::new();
    let handlers: Vec<Box<dyn Handler>> = Vec::new();
    let rng = SmallRng::seed_from_u64(1_234_567);
    (sys, x, swap, ws, handlers, rng)
}

fn phase_info() -> PhaseInfo {
    PhaseInfo {
        phase: 0,
        total_phases: 1,
        molecule_type: None,
    }
}

fn movebad_cfg() -> MoveBadConfig<'static> {
    MoveBadConfig {
        movefrac: 0.05,
        maxmove_per_type: &[],
        movebadrandom: false,
        gencan_maxit: 20,
    }
}

fn bench_run_iteration(c: &mut Criterion) {
    let mut group = c.benchmark_group("run_iteration");
    group.sample_size(10);
    group.measurement_time(Duration::from_millis(500));

    let pi = phase_info();
    let mb = movebad_cfg();
    let gp = GencanParams::default();

    group.bench_function("step", |b| {
        b.iter_batched(
            build_snapshot,
            |(mut sys, mut x, mut swap, mut ws, mut handlers, mut rng)| {
                let mut flast = 0.0_f64;
                let mut fimp_prev = F::INFINITY;
                let mut radscale = 1.0_f64;
                let out = run_iteration(
                    0,
                    10,
                    true,
                    0,
                    pi,
                    0.01,
                    true,
                    &mb,
                    &gp,
                    &mut sys,
                    &mut x,
                    &mut swap,
                    &mut flast,
                    &mut fimp_prev,
                    &mut radscale,
                    // Optimizer bindings are `pub(crate)`; an empty slice
                    // coerces without naming the type.
                    #[cfg(feature = "ff")]
                    &mut [],
                    &mut handlers,
                    &mut ws,
                    &mut rng,
                );
                debug_assert!(matches!(
                    out,
                    IterOutcome::Continue | IterOutcome::Converged | IterOutcome::EarlyStop
                ));
                std::hint::black_box(out);
            },
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_run_iteration);
criterion_main!(benches);
use rand::SeedableRng;
use rand::rngs::SmallRng;
