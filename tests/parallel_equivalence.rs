//! Numerical-equivalence test: `compute_fg` with
//! `PackContext::parallel_pair_eval = true` must produce the same `f`,
//! `fdist`, and full gradient as the serial path. Guards against any
//! float-summation reordering regression inside the rayon reduce +
//! parallel merge.
//!
//! Parallelism is user-opt-in via `Molpack::parallel_pair_eval(bool)`
//! → `PackContext::parallel_pair_eval`. These tests flip the flag
//! explicitly so the parallel code path is actually exercised instead
//! of relying on a size heuristic.

#![cfg(feature = "rayon")]

use molpack::objective::{compute_f, compute_fg};
use molpack::{F, PackContext};
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

/// Same synthetic water-in-a-box builder as `benches/pair_kernel.rs`.
/// Kept in-crate so the equivalence test uses the same workload shape
/// the bench already measures.
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
    sys.coor = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
    sys.radius.fill(1.0);
    sys.radius_ini.fill(1.0);
    sys.fscale.fill(1.0);

    for imol in 0..n_mols {
        for iatom in 0..atoms_per_mol {
            let icart = imol * atoms_per_mol + iatom;
            sys.ibtype[icart] = 0;
            sys.ibmol[icart] = imol;
        }
    }
    sys.iratom_offsets = vec![0; ntotat + 1];
    sys.iratom_data.clear();

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
    let _ = compute_f(&x, &mut sys);
    sys.reset_eval_counters();

    (sys, x)
}

/// `compute_fg` with `parallel_pair_eval = true` produces the same `f`
/// and gradient as the serial path (`compute_f` + `compute_g`). Both
/// consult the same inner pair kernel; differences can only come from
/// rayon's tree-order reduce in the fused path.
#[test]
fn compute_fg_parallel_matches_compute_g_serial_large_system() {
    let (mut sys_a, x) = build_water_box(1000, 80.0, 0xC0FFEE);
    let (mut sys_b, _) = build_water_box(1000, 80.0, 0xC0FFEE);
    // Explicit opt-in: exercise the parallel branch of
    // accumulate_pair_fg. `sys_b` stays serial (default false).
    sys_a.parallel_pair_eval = true;

    // Fused path (parallel — flag is on).
    let mut g_fused = vec![0.0; x.len()];
    let f_fused = compute_fg(&x, &mut sys_a, &mut g_fused);

    // Serial path via compute_f + compute_g; compute_g is not
    // parallelised, so its gradient is the reference ordering.
    let f_sep = compute_f(&x, &mut sys_b);
    let mut g_sep = vec![0.0; x.len()];
    molpack::objective::compute_g(&x, &mut sys_b, &mut g_sep);

    // f should be numerically identical up to the partial-sum
    // reordering inside rayon's reduce — both sum the same set of
    // scalar contributions, but parallel reduce aggregates them in
    // tree order. Allow 1e-8 relative tolerance.
    let rel_f = (f_fused - f_sep).abs() / f_sep.abs().max(1.0);
    assert!(
        rel_f < 1e-8,
        "f mismatch: fused={f_fused} sep={f_sep} rel={rel_f}"
    );

    // Gradient contributions for each atom pair are identical in both
    // paths; the parallel merge just adds them up in a different
    // chunk order. That is associative on f64 *modulo* rounding, so a
    // relative tolerance of 1e-9 is what to expect in practice.
    let mut max_abs_err: F = 0.0;
    let mut max_rel_err: F = 0.0;
    for (i, (&a, &b)) in g_fused.iter().zip(&g_sep).enumerate() {
        let abs = (a - b).abs();
        let rel = abs / a.abs().max(b.abs()).max(1.0);
        if abs > max_abs_err {
            max_abs_err = abs;
        }
        if rel > max_rel_err {
            max_rel_err = rel;
        }
        assert!(
            rel < 1e-9,
            "g[{i}] relative mismatch: fused={a} sep={b} rel={rel}"
        );
    }
    eprintln!("parallel_equivalence: max_abs={max_abs_err:e} max_rel={max_rel_err:e}");
}

/// Small-system path with parallelism still opted-in: even at a
/// workload size where parallel is net slower, the numerical result
/// must match the serial reference. Guards against regressions where
/// the parallel partition/merge subtly diverges (e.g. a stale partial
/// buffer slot) that a bigger workload might hide through noise.
#[test]
fn compute_fg_small_system_parallel_matches_serial() {
    let (mut sys_a, x) = build_water_box(50, 12.0, 0xBEEF);
    let (mut sys_b, _) = build_water_box(50, 12.0, 0xBEEF);
    sys_a.parallel_pair_eval = true;

    let mut g_fused = vec![0.0; x.len()];
    let f_fused = compute_fg(&x, &mut sys_a, &mut g_fused);

    let f_sep = compute_f(&x, &mut sys_b);
    let mut g_sep = vec![0.0; x.len()];
    molpack::objective::compute_g(&x, &mut sys_b, &mut g_sep);

    assert!((f_fused - f_sep).abs() < 1e-10);
    for (a, b) in g_fused.iter().zip(&g_sep) {
        assert!((a - b).abs() < 1e-10, "g mismatch: {a} vs {b}");
    }
}

// ── pack-level parity: full Molpack::pack run, serial vs parallel ──────────

/// The kernel-level tests above lock in `compute_fg` parity, but the
/// *full* pack driver layers gencan, relaxer, movebad, RNG sampling, and
/// handler IO on top. A subtle drift inside any of those that depends
/// on the parallel branch (e.g. an evaluation-order dependency in the
/// movebad heuristic) would not surface in `compute_fg` alone. This
/// test exercises a complete `Molpack::pack` run twice — same seed,
/// same targets — once serial and once with `parallel_pair_eval(true)`,
/// and asserts that the final atom coordinates agree to within
/// `1e-10`. Bit-exactness across rayon merge order is not guaranteed,
/// so we use an absolute tolerance rather than `==`.
///
/// Small workload (3 waters in a 30 Å box, 5 outer loops) keeps this
/// fast enough to live in the default tier (< 200 ms).
#[test]
fn pack_seed_parity_serial_vs_parallel() {
    use molpack::{InsideBoxRestraint, Molpack, Target};

    fn build_target() -> Target {
        let positions = vec![[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];
        let radii = vec![1.52, 1.20, 1.20];
        Target::from_coords(&positions, &radii, 3).with_restraint(InsideBoxRestraint::new(
            [0.0, 0.0, 0.0],
            [30.0, 30.0, 30.0],
            [false; 3],
        ))
    }

    const SEED: u64 = 0xA11CE;
    const MAX_LOOPS: usize = 5;

    let serial = Molpack::new()
        .with_seed(SEED)
        .with_parallel_eval(false)
        .pack(&[build_target()], MAX_LOOPS)
        .expect("serial pack failed");

    let parallel = Molpack::new()
        .with_seed(SEED)
        .with_parallel_eval(true)
        .pack(&[build_target()], MAX_LOOPS)
        .expect("parallel pack failed");

    assert_eq!(
        serial.natoms(),
        parallel.natoms(),
        "atom count diverged: {} serial vs {} parallel",
        serial.natoms(),
        parallel.natoms()
    );

    let s = serial.positions();
    let p = parallel.positions();
    let mut max_err: F = 0.0;
    for (i, (a, b)) in s.iter().zip(p.iter()).enumerate() {
        for k in 0..3 {
            let err = (a[k] - b[k]).abs();
            if err > max_err {
                max_err = err;
            }
            assert!(
                err < 1e-10,
                "atom {i} axis {k}: serial={} parallel={} err={err}",
                a[k],
                b[k]
            );
        }
    }
    eprintln!("pack_seed_parity: max coord diff = {max_err:e}");
}
