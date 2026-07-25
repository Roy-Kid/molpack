//! Regression microbench for the fused pair/objective kernel.
//!
//! Times `objective::compute_f` (function only) and `objective::compute_fg`
//! (the fused function+gradient pass GENCAN drives on every inner iteration)
//! on one small, self-synthesized water-box snapshot — the same state the
//! outer GENCAN loop hands the objective on every evaluate.
//!
//! This is a **regression alarm**, not a throughput measurement: it uses one
//! small workload (30 water-like molecules = 90 atoms) so it runs in well
//! under a second and catches a kernel regression, not peak scaling. The
//! larger, size-swept scaling study lives in `examples/mt_scaling`.
//!
//! Geometry is synthesized in-process (no external test data) so the bench
//! needs no `io` feature.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use molpack::objective::{compute_f, compute_fg};
use molpack::{F, PackContext};

/// Small representative size: keeps the whole bench << 1 s.
const N_MOLS: usize = 30;

/// Build a synthetic `PackContext` of `n_mols` water-like molecules
/// (atoms_per_mol = 3) jittered inside a cubic box, no PBC wrap, no short
/// radius, no fixed atoms. Returns the context plus a packer-format `x`
/// vector (`[COM₀..COMₙ, euler₀..eulerₙ]`). Mirrors the synthesis in
/// `examples/mt_scaling/main.rs`.
fn build_water_box(n_mols: usize, box_side: F, seed: u64) -> (PackContext, Vec<F>) {
    let atoms_per_mol = 3usize;
    let ntotat = n_mols * atoms_per_mol;
    let ntype = 1usize;
    let mut sys = PackContext::new(ntotat, n_mols, ntype);
    sys.ntype_with_fixed = ntype;

    sys.nmols = vec![n_mols];
    sys.natoms = vec![atoms_per_mol];
    sys.idfirst = vec![0];
    sys.comptype = vec![true; ntype];
    sys.constrain_rot = vec![[false; 3]; ntype];
    sys.rot_bound = vec![[[0.0; 2]; 3]; ntype];

    // Reference coords: water-like O + 2H layout.
    sys.coor = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

    // Radii: tolerance/2 = 1.0 Å for tolerance = 2.0.
    sys.radius.fill(1.0);
    sys.radius_ini.fill(1.0);
    sys.fscale.fill(1.0);

    // Per-atom molecule/type tags (same layout as `Molpack::pack` writes).
    for imol in 0..n_mols {
        for iatom in 0..atoms_per_mol {
            let icart = imol * atoms_per_mol + iatom;
            sys.ibtype[icart] = 0;
            sys.ibmol[icart] = imol;
        }
    }

    // Restraint CSR — no restraints, so all offsets = 0.
    sys.iratom_offsets = vec![0; ntotat + 1];
    sys.iratom_data.clear();

    // Cell geometry over a padded, non-periodic box.
    let pad: F = 3.0;
    let side = box_side + 2.0 * pad;
    let origin = [-pad, -pad, -pad];
    let mut origin_arr = molrs::types::F3::zeros(3);
    for k in 0..3 {
        origin_arr[k] = origin[k];
    }
    sys.simbox = molrs::spatial::region::simbox::SimBox::cube(side, origin_arr, [false; 3])
        .expect("bench cell");
    let cell_side: F = 2.0; // ≈ 1.01 * 2*radius_ini
    sys.grid = molrs::spatial::neighbors::CellGrid::for_cutoff(&sys.simbox, cell_side);
    sys.resize_cell_arrays();

    sys.sizemin = origin;
    sys.sizemax = [origin[0] + side, origin[1] + side, origin[2] + side];

    sys.sync_atom_props();

    // Random COM + Euler angles (x layout = [COM₀..COMₙ, euler₀..eulerₙ]).
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut x = vec![0.0 as F; 6 * n_mols];
    for imol in 0..n_mols {
        x[3 * imol] = rng.random::<F>() * box_side;
        x[3 * imol + 1] = rng.random::<F>() * box_side;
        x[3 * imol + 2] = rng.random::<F>() * box_side;
        let base = 3 * n_mols + 3 * imol;
        x[base] = rng.random::<F>() * std::f64::consts::TAU;
        x[base + 1] = rng.random::<F>() * std::f64::consts::TAU;
        x[base + 2] = rng.random::<F>() * std::f64::consts::TAU;
    }

    // Warm up: run compute_f once to populate cells + active_cells.
    let _ = compute_f(&x, &mut sys);
    sys.reset_eval_counters();

    (sys, x)
}

fn bench_pair_kernel(c: &mut Criterion) {
    let mut group = c.benchmark_group("pair_kernel");
    group.sample_size(10);
    group.measurement_time(Duration::from_millis(500));

    let box_side = 30.0 + (N_MOLS as F).cbrt() * 2.5;

    // compute_f: constant `x`, so every call repeats the same work.
    {
        let (mut sys, x) = build_water_box(N_MOLS, box_side, 0x5eed);
        group.bench_function(format!("compute_f/atoms={}", N_MOLS * 3), |b| {
            b.iter(|| {
                std::hint::black_box(compute_f(&x, &mut sys));
            });
        });
    }

    // compute_fg: fused function + gradient (the GENCAN inner-loop kernel).
    {
        let (mut sys, x) = build_water_box(N_MOLS, box_side, 0x5eed);
        let mut g = vec![0.0; x.len()];
        group.bench_function(format!("compute_fg/atoms={}", N_MOLS * 3), |b| {
            b.iter(|| {
                std::hint::black_box(compute_fg(&x, &mut sys, &mut g));
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_pair_kernel);
criterion_main!(benches);
use rand::RngExt;
use rand::SeedableRng;
use rand::rngs::SmallRng;
