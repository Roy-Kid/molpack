//! Hot-loop microbench for the pair-evaluation kernel.
//!
//! Measures `objective::compute_f` and `objective::compute_fg` on synthetic
//! mid-pack snapshots of varying size. Each iteration calls `compute_f` /
//! `compute_fg` on a prebuilt [`PackContext`] whose cells are already sized
//! and xcart positions already jittered — i.e. the same state the outer
//! GENCAN loop hands to the objective on every evaluate.
//!
//! Unlike `objective_dispatch` / `run_iteration` (which use `ntotat=4` to
//! measure function-call edge costs), this bench uses system sizes where
//! the per-pair kernel actually dominates (≥ 3000 atoms) so that AoS
//! packing, PBC short-circuit, and branch-elimination wins show up.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use molpack::objective::{compute_f, compute_fg};
use molpack::{F, PackContext};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

/// Build a synthetic PackContext representative of the pair-kernel hot path:
/// `n_mols` water-like molecules (atoms_per_mol=3) jittered inside a cubic
/// box of side `box_side` Å, no PBC wrap, no short radius, no fixed atoms.
/// Returns the context plus an `x` vector (3N COM + 3N Euler) that is the
/// packer's canonical input format.
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

    // Radii: Packmol tolerance/2 = 1.0 Å for tolerance=2.0.
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

    // Cell geometry.  We use a padded box so `setcell` does no wrapping and
    // the PBC branch short-circuits (matches `Molpack::pack` without `.pbc()`).
    let pad: F = 3.0;
    sys.pbc_min = [-pad, -pad, -pad];
    sys.pbc_length = [box_side + 2.0 * pad; 3];
    let cell_side: F = 2.0; // ≈ 1.01 * 2*radius_ini
    for k in 0..3 {
        sys.ncells[k] = ((sys.pbc_length[k] / cell_side).floor() as usize).max(1);
        sys.cell_length[k] = sys.pbc_length[k] / sys.ncells[k] as F;
    }
    sys.resize_cell_arrays();

    sys.sizemin = sys.pbc_min;
    sys.sizemax = [
        sys.pbc_min[0] + sys.pbc_length[0],
        sys.pbc_min[1] + sys.pbc_length[1],
        sys.pbc_min[2] + sys.pbc_length[2],
    ];

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

fn bench_compute_f(c: &mut Criterion) {
    let mut group = c.benchmark_group("pair_kernel/compute_f");
    group.sample_size(30);
    for &n in &[200usize, 1000, 3000] {
        let label = format!("n={}_atoms={}", n, n * 3);
        let (sys_template, x) = build_water_box(n, 30.0 + (n as F).cbrt() * 2.5, 0x5eed);
        group.bench_function(&label, |b| {
            b.iter_batched_ref(
                || (sys_template_clone(&sys_template), x.clone()),
                |(sys, x)| {
                    std::hint::black_box(compute_f(x, sys));
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

fn bench_compute_fg(c: &mut Criterion) {
    let mut group = c.benchmark_group("pair_kernel/compute_fg");
    group.sample_size(30);
    for &n in &[200usize, 1000, 3000] {
        let label = format!("n={}_atoms={}", n, n * 3);
        let (sys_template, x) = build_water_box(n, 30.0 + (n as F).cbrt() * 2.5, 0x5eed);
        let g_len = 6 * n;
        group.bench_function(&label, |b| {
            b.iter_batched_ref(
                || (sys_template_clone(&sys_template), x.clone(), vec![0.0; g_len]),
                |(sys, x, g)| {
                    std::hint::black_box(compute_fg(x, sys, g));
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

/// Custom clone for PackContext — only the fields the pair-kernel bench
/// reads or mutates.  Full `Clone` isn't derived because a handful of
/// internal types (e.g. `Frame`, `Constraints`) are intentionally
/// non-`Clone`; we fake it by rebuilding just the subset we need.
fn sys_template_clone(src: &PackContext) -> PackContext {
    let ntotat = src.ntotat;
    let ntotmol = src.ntotmol;
    let ntype = src.ntype;
    let mut dst = PackContext::new(ntotat, ntotmol, ntype);
    dst.ntype_with_fixed = src.ntype_with_fixed;
    dst.nmols = src.nmols.clone();
    dst.natoms = src.natoms.clone();
    dst.idfirst = src.idfirst.clone();
    dst.comptype = src.comptype.clone();
    dst.constrain_rot = src.constrain_rot.clone();
    dst.rot_bound = src.rot_bound.clone();
    dst.coor = src.coor.clone();
    dst.radius = src.radius.clone();
    dst.radius_ini = src.radius_ini.clone();
    dst.fscale = src.fscale.clone();
    dst.ibtype = src.ibtype.clone();
    dst.ibmol = src.ibmol.clone();
    dst.fixedatom = src.fixedatom.clone();
    dst.iratom_offsets = src.iratom_offsets.clone();
    dst.iratom_data = src.iratom_data.clone();
    dst.pbc_min = src.pbc_min;
    dst.pbc_length = src.pbc_length;
    dst.ncells = src.ncells;
    dst.cell_length = src.cell_length;
    dst.sizemin = src.sizemin;
    dst.sizemax = src.sizemax;
    dst.resize_cell_arrays();
    dst.sync_atom_props();
    dst
}

criterion_group!(benches, bench_compute_f, bench_compute_fg);
criterion_main!(benches);
