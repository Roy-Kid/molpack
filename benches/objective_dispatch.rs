//! Regression microbench for `Objective::evaluate` dispatch cost.
//!
//! A pack run drives a `&mut dyn Objective` through GENCAN; every
//! `sys.evaluate(...)` call site pays one vtable indirection. This bench
//! isolates that dispatch cost: the same trivial evaluation through the
//! inherent `PackContext::evaluate` impl (`via_inherent`) versus through a
//! `&mut dyn Objective` receiver (`via_dyn`).
//!
//! Setup mirrors the old A.6 dispatch bench: an empty-molecule `PackContext`
//! (`ntotat = 4`, `ntotmol = 0`), `FOnly` mode, empty gradient slot — so
//! `evaluate` returns `(0, 0, 0)` and the bench measures boundary cost, not
//! body cost. Tiny and self-contained: no external test data, no features.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use molpack::{EvalMode, F, Objective, PackContext};

fn build_ctx() -> PackContext {
    let ntotat = 4;
    let mut sys = PackContext::new(ntotat, 0, 0);
    sys.radius.fill(0.75);
    sys.radius_ini.fill(1.5);
    sys.work.radiuswork.resize(ntotat, 0.0);
    sys.sync_atom_props();
    sys
}

fn bench_objective_dispatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("objective_dispatch");
    group.sample_size(10);
    group.measurement_time(Duration::from_millis(500));

    let x: Vec<F> = Vec::new();

    let mut sys = build_ctx();
    group.bench_function("via_inherent", |b| {
        b.iter(|| {
            let out = PackContext::evaluate(&mut sys, &x, EvalMode::FOnly, None);
            std::hint::black_box(out);
        });
    });

    let mut sys = build_ctx();
    group.bench_function("via_dyn", |b| {
        b.iter(|| {
            let obj: &mut dyn Objective = &mut sys;
            let out = obj.evaluate(&x, EvalMode::FOnly, None);
            std::hint::black_box(out);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_objective_dispatch);
criterion_main!(benches);
