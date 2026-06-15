//! Parallel-scaling harness for the pair-evaluation kernel on a **mixture**
//! workload (two species: an 8-atom rigid molecule + a 3-atom water-like
//! molecule, 1:2 count ratio — the same composition as the molpack-jcc
//! `bench_mixture.py` paper bench).
//!
//! For each system size it times one `compute_fg` evaluation (the fused
//! function+gradient pass GENCAN drives on every inner iteration) using the
//! *same* rayon kernel at 1/2/4/8 worker threads — the 1-worker run is the
//! serial baseline, so the speed-up isolates thread scaling rather than mixing
//! in the separate (phase-structured) serial code path — and prints
//!
//! ```text
//! atoms,kernel,mode,threads,us_per_eval,speedup_vs_serial
//! ```
//!
//! to stdout — the schema `molpack-jcc/bench/fig_mt_scaling.py` reads for
//! `fig_parallel_kernel.{pdf,png}`. The default size ladder reaches the
//! million-atom regime, where the parallel kernel dominates and scaling
//! approaches the machine ceiling.
//!
//! Run (needs the rayon feature):
//! ```sh
//! cargo run --release --example mt_scaling --features rayon
//! cargo run --release --example mt_scaling --features rayon -- 60000 1000000
//! ```
//!
//! Each `compute_fg` is fed a freshly jittered `x`, so the geometry cache
//! misses every call and the full expand→pair→project path runs — matching
//! what GENCAN's line search actually exercises.

use std::time::{Duration, Instant};

use molpack::objective::compute_fg;
use molpack::{F, PackContext};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

/// Number density (atoms / Å³) held constant across sizes so the comparison is
/// weak-scaling: cell occupancy — and therefore work-per-cell — stays fixed and
/// only the cell *count* grows with N. ~0.06 is a representative mid-pack
/// density for a tolerance-2.0 system (radius 1.0 Å).
const NUMBER_DENSITY: F = 0.06;
/// rayon thread counts to sweep.
const THREADS: [usize; 4] = [1, 2, 4, 8];
/// Wall-time spent in each timing repetition.
const TIME_BUDGET: Duration = Duration::from_millis(150);
/// Repetitions per (size, threads) point. Each rep times the *same* kernel at
/// every worker count back-to-back and forms a within-rep speed-up ratio against
/// the single-worker time; the reported number is the **median** of those ratios
/// across reps — robust to both the throttled-baseline outliers that inflate a
/// `max` and the throttled-parallel outliers that deflate a `min`.
const REPS: usize = 7;
/// Busy warm-up run (untimed) before timing each config. It keeps the worker
/// thread(s) resident on performance cores and ramps the clock up — crucially,
/// it avoids the efficiency-core demotion an *idle* pause would trigger on the
/// asymmetric Apple-silicon scheduler, which otherwise inflates the
/// single-thread serial baseline. No shell/idle sleep is used anywhere.
const WARMUP: Duration = Duration::from_millis(80);

/// Reference coords for the 8-atom rigid species (centered, ~5 Å across).
fn mol_a_coords() -> Vec<[F; 3]> {
    let mut rng = SmallRng::seed_from_u64(0xA);
    let mut c: Vec<[F; 3]> = (0..8)
        .map(|_| {
            [
                rng.random::<F>() * 2.0 - 1.0,
                rng.random::<F>() * 2.0 - 1.0,
                rng.random::<F>() * 2.0 - 1.0,
            ]
        })
        .collect();
    // Center at origin.
    let mut mean = [0.0; 3];
    for p in &c {
        for k in 0..3 {
            mean[k] += p[k] / c.len() as F;
        }
    }
    for p in &mut c {
        for k in 0..3 {
            p[k] = (p[k] - mean[k]) * 2.5;
        }
    }
    c
}

/// Reference coords for the 3-atom water-like species.
fn mol_b_coords() -> Vec<[F; 3]> {
    vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]
}

/// Build a two-species mixture `PackContext` of ~`target_atoms` atoms, sized to
/// hold them at [`NUMBER_DENSITY`]. `n_b = 2 * n_a` (1:2 count ratio).
/// Returns the context plus a packer-format `x` (`[COM₀..; euler₀..]`).
fn build_mixture(target_atoms: usize, seed: u64) -> (PackContext, Vec<F>) {
    let a = mol_a_coords();
    let b = mol_b_coords();
    let (na_at, nb_at) = (a.len(), b.len()); // 8, 3
    let per_unit = na_at + 2 * nb_at; // atoms per (1 A + 2 B) unit = 14
    let units = (target_atoms / per_unit).max(1);
    let (n_a, n_b) = (units, 2 * units);
    let ntotmol = n_a + n_b;
    let ntotat = n_a * na_at + n_b * nb_at;
    let ntype = 2usize;

    let box_side = (ntotat as F / NUMBER_DENSITY).cbrt();

    let mut sys = PackContext::new(ntotat, ntotmol, ntype);
    sys.ntype_with_fixed = ntype;
    sys.nmols = vec![n_a, n_b];
    sys.natoms = vec![na_at, nb_at];
    sys.idfirst = vec![0, na_at];
    sys.comptype = vec![true; ntype];
    sys.constrain_rot = vec![[false; 3]; ntype];
    sys.rot_bound = vec![[[0.0; 2]; 3]; ntype];

    let mut coor = a.clone();
    coor.extend_from_slice(&b);
    sys.coor = coor;

    // Radii: tolerance/2 = 1.0 for tolerance 2.0, as Molpack::pack assigns.
    sys.radius.fill(1.0);
    sys.radius_ini.fill(1.0);
    sys.fscale.fill(1.0);

    // Per-atom (itype, imol-within-type) tags, in the packer's atom layout:
    // all of type 0 first, then all of type 1.
    let mut icart = 0usize;
    for imol in 0..n_a {
        for _ in 0..na_at {
            sys.ibtype[icart] = 0;
            sys.ibmol[icart] = imol;
            icart += 1;
        }
    }
    for imol in 0..n_b {
        for _ in 0..nb_at {
            sys.ibtype[icart] = 1;
            sys.ibmol[icart] = imol;
            icart += 1;
        }
    }

    // No restraints: empty CSR.
    sys.iratom_offsets = vec![0; ntotat + 1];
    sys.iratom_data.clear();

    // Cell grid over a padded box (so setcell never wraps; non-PBC path).
    let pad: F = 3.0;
    sys.pbc_min = [-pad, -pad, -pad];
    sys.pbc_length = [box_side + 2.0 * pad; 3];
    let cell_side: F = 2.0;
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

    // Random COM + Euler (x layout = [COM for all mols; euler for all mols]).
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut x = vec![0.0 as F; 6 * ntotmol];
    for imol in 0..ntotmol {
        x[3 * imol] = rng.random::<F>() * box_side;
        x[3 * imol + 1] = rng.random::<F>() * box_side;
        x[3 * imol + 2] = rng.random::<F>() * box_side;
        let base = 3 * ntotmol + 3 * imol;
        x[base] = rng.random::<F>() * std::f64::consts::TAU;
        x[base + 1] = rng.random::<F>() * std::f64::consts::TAU;
        x[base + 2] = rng.random::<F>() * std::f64::consts::TAU;
    }

    // Warm up the cells / active_cells.
    let mut g = vec![0.0; x.len()];
    let _ = compute_fg(&x, &mut sys, &mut g);
    sys.reset_eval_counters();
    (sys, x)
}

/// One timing repetition: run `compute_fg` (the fused function+gradient pass
/// GENCAN drives on every inner iteration) until [`TIME_BUDGET`] elapses,
/// jittering one coordinate each call so the geometry cache misses and the full
/// kernel runs every time. Returns mean microseconds per evaluation for this
/// rep. The budget is checked every iteration (not every Nth) so a single
/// million-atom evaluation, which already exceeds the budget on its own, does
/// not force a long minimum batch.
fn time_rep(sys: &mut PackContext, x: &mut [F], g: &mut [F]) -> F {
    let mut iters = 0u64;
    // jitter step chosen tiny so the system barely drifts over the run.
    let step = 1e-7;
    let t0 = Instant::now();
    loop {
        x[0] += step;
        std::hint::black_box(compute_fg(x, sys, g));
        iters += 1;
        if t0.elapsed() >= TIME_BUDGET {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    (secs / iters as f64) * 1e6
}

/// Per size: time the *same* kernel at every worker count **back-to-back within
/// each rep** so they share the same thermal state, then report, per worker
/// count, the *median* within-rep speed-up against the single-worker time across
/// [`REPS`] reps. Taking the ratio from the same rep cancels the sustained-load
/// throttling that otherwise depresses many-thread numbers (one worker heats the
/// chip far less than an 8-worker run), making the scaling curve reproducible
/// across invocations.
///
/// The baseline is the single-worker run of the **identical** kernel, not the
/// separate serial code path, so the ratio reflects parallelism only — an
/// algorithm change can never leak into the reported speed-up.
///
/// A **single** context is shared across all configs (only the rayon pool
/// differs) — at million-atom sizes one context's cell lists and gradient
/// buffers already run to ~1 GB, so building one per worker count would not fit.
/// The continuously jittered `x` drifts by < 1e-3 Å over a run, negligible for
/// the kernel cost.
fn run_size(target_atoms: usize) {
    let (mut sys, mut x) = build_mixture(target_atoms, 0x5eed);
    let ntotat = sys.ntotat;

    let pools: Vec<rayon::ThreadPool> = THREADS
        .iter()
        .map(|&nt| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(nt)
                .build()
                .expect("build rayon pool")
        })
        .collect();

    let mut base_us: Vec<F> = Vec::with_capacity(REPS);
    let mut speedups: Vec<Vec<F>> = vec![Vec::with_capacity(REPS); THREADS.len()];
    let mut par_us: Vec<Vec<F>> = vec![Vec::with_capacity(REPS); THREADS.len()];

    // CRITICAL FAIRNESS: every point — including the serial baseline — runs the
    // *identical* kernel (`parallel_pair_eval = true`). Only the rayon worker
    // count differs, so the ratio measures thread scaling alone, never an
    // algorithm change. (Using the serial code path `parallel_pair_eval = false`
    // as the baseline would compare two *different* implementations: the
    // single-worker run would then be spuriously faster/slower than serial for
    // reasons unrelated to parallelism.) Warm up each pool size so its workers
    // are P-core-resident and clocked up before the first timed rep.
    sys.parallel_pair_eval = true;
    for pool in &pools {
        pool.install(|| warmup(&mut sys, &mut x));
    }

    for _ in 0..REPS {
        // Time the same kernel at every worker count back-to-back so all points
        // share the same thermal state. THREADS[0] == 1 is the single-worker
        // baseline: same code path, no parallelism, so its self-ratio is 1.0 by
        // construction (strong-scaling convention S(1) = 1).
        sys.parallel_pair_eval = true;
        let mut t: Vec<F> = Vec::with_capacity(pools.len());
        for pool in &pools {
            t.push(pool.install(|| time_kernel_one(&mut sys, &mut x)));
        }
        let base = t[0];
        base_us.push(base);
        for (i, &ti) in t.iter().enumerate() {
            speedups[i].push(base / ti); // within-rep ratio → thermally fair
            par_us[i].push(ti);
        }
    }

    // CSV schema kept stable (`kernel` column = "fg") for the figure script. The
    // `serial` row is the single-worker run of the *same* parallel kernel, so
    // every speed-up isolates thread scaling.
    println!("{ntotat},fg,serial,1,{:.1},1.00", median(&mut base_us));
    for (i, &nt) in THREADS.iter().enumerate() {
        println!(
            "{ntotat},fg,parallel,{nt},{:.1},{:.2}",
            median(&mut par_us[i]),
            median(&mut speedups[i])
        );
    }
}

/// Median of a slice (sorts in place). Average of the two middle elements for
/// even lengths.
fn median(v: &mut [F]) -> F {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// One timing rep (see [`time_rep`]).
fn time_kernel_one(sys: &mut PackContext, x: &mut [F]) -> F {
    let mut g = vec![0.0; x.len()];
    time_rep(sys, x, &mut g)
}

/// Busy warm-up: drive `compute_fg` for [`WARMUP`] without recording timings.
fn warmup(sys: &mut PackContext, x: &mut [F]) {
    let mut g = vec![0.0; x.len()];
    let step = 1e-7;
    let t0 = Instant::now();
    loop {
        x[0] += step;
        std::hint::black_box(compute_fg(x, sys, &mut g));
        if t0.elapsed() >= WARMUP {
            break;
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    // Default ladder is log-spaced from a few hundred atoms to a million, so the
    // speed-up-vs-size curve shows both the small-system regime where the rayon
    // overhead makes parallel *slower* than serial and the large-system plateau
    // where the speed-up saturates.
    let sizes: Vec<usize> = if args.is_empty() {
        vec![
            300, 1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000,
        ]
    } else {
        args.iter()
            .map(|s| s.parse().expect("size must be an integer"))
            .collect()
    };

    println!("atoms,kernel,mode,threads,us_per_eval,speedup_vs_serial");
    for &s in &sizes {
        run_size(s);
    }
}
