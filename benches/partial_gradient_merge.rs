//! Microbench for the parallel-pair-kernel's gradient merge step.
//!
//! `accumulate_pair_fg_parallel` accumulates per-thread partial
//! gradient buffers then reduces them into `sys.work.gxcar`. Before
//! Fix 3 that reduction was a serial nested loop over
//! `n_threads × ntotat × 3` scalar adds — O(n_threads) factor Amdahl'd
//! against the parallel pair work. This bench times just the merge step
//! with `merge_partials_serial` (the pre-fix version) and
//! `merge_partials_parallel` (the post-fix version) across the `ntotat`
//! range where the pack kernel actually operates: from ~1k (small
//! mixtures) up to 64k atoms (`pack_spherical`-class).
//!
//! Used to pick `MERGE_PARALLEL_THRESHOLD_NTOTAT` in
//! `work_buffers.rs`. Run with:
//!
//! ```sh
//! cargo bench --bench partial_gradient_merge --features rayon
//! ```

#![cfg(feature = "rayon")]

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use molpack::F;
use molpack::context::WorkBuffers;

/// Build `n_threads` partial buffers of `ntotat` triples, prefilled
/// with a deterministic non-zero pattern so the merge actually has
/// work (and so compilers can't elide the adds as dead code).
fn build_partials(n_threads: usize, ntotat: usize) -> Vec<Vec<[F; 3]>> {
    (0..n_threads)
        .map(|t| {
            (0..ntotat)
                .map(|i| {
                    let base = (t * 7 + i) as F;
                    [base * 1e-3, base * 2e-3, base * 3e-3]
                })
                .collect()
        })
        .collect()
}

fn bench_merge(c: &mut Criterion) {
    let mut group = c.benchmark_group("partial_gradient_merge");
    group.sample_size(30);
    let n_threads = rayon::current_num_threads().max(1);

    for &ntotat in &[1_024usize, 4_096, 16_384, 65_536] {
        let partials = build_partials(n_threads, ntotat);

        group.bench_function(format!("serial/ntotat={}", ntotat), |b| {
            b.iter_batched_ref(
                || vec![[0.0 as F; 3]; ntotat],
                |out| {
                    WorkBuffers::merge_partials_serial(&partials, out);
                    std::hint::black_box(&out[0]);
                },
                BatchSize::SmallInput,
            );
        });

        group.bench_function(format!("parallel/ntotat={}", ntotat), |b| {
            b.iter_batched_ref(
                || vec![[0.0 as F; 3]; ntotat],
                |out| {
                    WorkBuffers::merge_partials_parallel(&partials, out);
                    std::hint::black_box(&out[0]);
                },
                BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, bench_merge);
criterion_main!(benches);
