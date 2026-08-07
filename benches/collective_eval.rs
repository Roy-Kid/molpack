//! Cost of one collective-restraint gradient evaluation as the species grows.
//!
//! `accumulate_collective_fg` runs on every GENCAN gradient evaluation, so its
//! per-evaluation cost is multiplied by thousands. It used to stage the group's
//! coordinates and gradients through two fresh `Vec`s; this bench is what makes
//! that visible.

use std::sync::Arc;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use molpack::F;
use molpack::context::PackContext;
use molpack::objective::compute_g;
use molpack::restraint::GaussianPlane;

fn system(nmol: usize) -> (PackContext, Vec<F>) {
    let mut sys = PackContext::new(nmol, nmol, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![nmol];
    sys.natoms = vec![1];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0]];
    sys.radius = vec![1.0; nmol];
    sys.radius_ini = vec![1.0; nmol];
    sys.fscale = vec![1.0; nmol];
    for i in 0..nmol {
        sys.ibmol[i] = i;
    }
    sys.sync_atom_props();
    sys.init1 = false;

    // A real cell partition: without one every atom lands in a single cell and
    // the pair loop degrades to O(N^2), which would swamp what this bench is
    // about.
    let side = 12.0 * (nmol as F).cbrt();
    sys.simbox = molrs::spatial::simbox::SimBox::cube(side, molrs::types::F3::zeros(3), [false; 3])
        .expect("cell");
    sys.grid = molrs::spatial::neighbors::CellGrid::for_cutoff(&sys.simbox, 2.5);
    sys.resize_cell_arrays();
    sys.collective = vec![(
        0usize,
        Arc::new(GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1.0, 0.0, 20.0)) as Arc<_>,
    )];

    // Spread the molecules through the cell so the pair loop is representative.
    let mut x = vec![0.0; 6 * nmol];
    let per_axis = (nmol as F).cbrt().ceil() as usize;
    for m in 0..nmol {
        let (i, j, k) = (
            m % per_axis,
            (m / per_axis) % per_axis,
            m / (per_axis * per_axis),
        );
        x[3 * m] = i as F * 12.0;
        x[3 * m + 1] = j as F * 12.0;
        x[3 * m + 2] = k as F * 12.0;
    }
    (sys, x)
}

fn bench_collective(c: &mut Criterion) {
    let mut group = c.benchmark_group("collective_eval");
    group.sample_size(10);
    group.measurement_time(Duration::from_millis(500));

    for &n in &[1_000usize, 10_000] {
        let (mut sys, x) = system(n);
        let mut g = vec![0.0; x.len()];
        group.bench_with_input(BenchmarkId::new("compute_g", n), &n, |b, _| {
            b.iter(|| {
                compute_g(&x, &mut sys, &mut g);
                std::hint::black_box(&g);
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_collective);
criterion_main!(benches);
