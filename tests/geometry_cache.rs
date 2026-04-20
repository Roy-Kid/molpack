#![allow(clippy::needless_range_loop)]
//! Tests for the geometry cache fast path in `compute_f` / `compute_fg` /
//! `compute_g`. The cache is hit when the caller evaluates at the same `x`
//! (and identical comptype / cell grid) as the previous call — in that case
//! the Cartesian expansion and cell-list rebuild are skipped and only the
//! pair / constraint kernels re-run on the stored state.
//!
//! These tests assert the cache path produces bit-identical results to the
//! fresh-rebuild path across several call sequences used by the packer.

use std::sync::Arc;

use molpack::objective::{compute_f, compute_fg};
use molpack::{F, InsideBoxRestraint, PackContext};

// ── setup helpers (mirror tests/gradient.rs patterns) ──────────────────────

fn setup_cells(sys: &mut PackContext, cell_n: usize, cell_len: F) {
    sys.ncells = [cell_n, cell_n, cell_n];
    sys.cell_length = [cell_len; 3];
    sys.pbc_min = [0.0; 3];
    sys.pbc_length = [cell_len * cell_n as F; 3];
    sys.resize_cell_arrays();
}

/// Three single-atom molecules inside a 5³ box with a pair-overlap setup.
fn mixed_system() -> (PackContext, Vec<F>) {
    let mut sys = PackContext::new(3, 3, 1);
    sys.ntype_with_fixed = 1;
    sys.nmols = vec![3];
    sys.natoms = vec![1];
    sys.idfirst = vec![0];
    sys.comptype = vec![true];
    sys.coor = vec![[0.0, 0.0, 0.0]];

    sys.radius = vec![1.0; 3];
    sys.radius_ini = vec![1.0; 3];
    sys.fscale = vec![1.0; 3];
    sys.ibmol = vec![0, 1, 2];
    sys.sync_atom_props();

    sys.restraints = vec![Arc::new(InsideBoxRestraint::new(
        [0.0, 0.0, 0.0],
        [5.0, 5.0, 5.0],
        [false; 3],
    ))];
    sys.iratom_offsets = vec![0, 1, 2, 3];
    sys.iratom_data = vec![0, 0, 0];

    setup_cells(&mut sys, 1, 10.0);

    // x = [com0(3), com1(3), com2(3), euler0(3), euler1(3), euler2(3)]
    let x = vec![
        6.0, 2.0, 2.0, // com0: outside box on +x, forces restraint penalty
        3.0, 2.0, 2.0, // com1: close to com2, forces pair penalty
        3.5, 2.5, 2.0, // com2
        0.1, 0.2, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, 0.4,
    ];
    (sys, x)
}

fn force_cache_miss(sys: &mut PackContext) {
    sys.work.cached_geometry_valid = false;
}

// ── compute_f cache ────────────────────────────────────────────────────────

#[test]
fn compute_f_cached_path_matches_fresh() {
    let (mut sys_a, x) = mixed_system();
    let (mut sys_b, _) = mixed_system();

    // Prime both systems: one cache, one fresh each call.
    let f_a1 = compute_f(&x, &mut sys_a);
    let fdist_a1 = sys_a.fdist;
    let frest_a1 = sys_a.frest;

    force_cache_miss(&mut sys_b);
    let f_b1 = compute_f(&x, &mut sys_b);
    let fdist_b1 = sys_b.fdist;
    let frest_b1 = sys_b.frest;

    assert_eq!(
        f_a1.to_bits(),
        f_b1.to_bits(),
        "compute_f first call must agree bitwise"
    );
    assert_eq!(fdist_a1.to_bits(), fdist_b1.to_bits());
    assert_eq!(frest_a1.to_bits(), frest_b1.to_bits());

    // Second call at same x: a hits cache, b forced miss.
    let f_a2 = compute_f(&x, &mut sys_a);
    force_cache_miss(&mut sys_b);
    let f_b2 = compute_f(&x, &mut sys_b);

    assert_eq!(f_a2.to_bits(), f_b2.to_bits());
    assert_eq!(f_a2.to_bits(), f_a1.to_bits(), "cache must be pure");
    assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
    assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
}

// ── compute_fg cache ───────────────────────────────────────────────────────

#[test]
fn compute_fg_cached_path_matches_fresh() {
    let (mut sys_a, x) = mixed_system();
    let (mut sys_b, _) = mixed_system();

    let mut g_a1 = vec![0.0; x.len()];
    let mut g_b1 = vec![0.0; x.len()];
    let f_a1 = compute_fg(&x, &mut sys_a, &mut g_a1);
    force_cache_miss(&mut sys_b);
    let f_b1 = compute_fg(&x, &mut sys_b, &mut g_b1);

    assert_eq!(f_a1.to_bits(), f_b1.to_bits());
    for i in 0..x.len() {
        assert_eq!(
            g_a1[i].to_bits(),
            g_b1[i].to_bits(),
            "compute_fg first call: g[{i}] mismatch {} vs {}",
            g_a1[i],
            g_b1[i]
        );
    }

    // Second call at same x — a hits cache, b forced miss.
    let mut g_a2 = vec![0.0; x.len()];
    let mut g_b2 = vec![0.0; x.len()];
    let f_a2 = compute_fg(&x, &mut sys_a, &mut g_a2);
    force_cache_miss(&mut sys_b);
    let f_b2 = compute_fg(&x, &mut sys_b, &mut g_b2);

    assert_eq!(f_a2.to_bits(), f_b2.to_bits());
    assert_eq!(f_a2.to_bits(), f_a1.to_bits(), "cache must be pure");
    for i in 0..x.len() {
        assert_eq!(g_a2[i].to_bits(), g_b2[i].to_bits());
        assert_eq!(g_a1[i].to_bits(), g_a2[i].to_bits());
    }
}

// ── cross-mode cache reuse (compute_fg → compute_f at same x) ──────────────

#[test]
fn compute_f_reuses_compute_fg_geometry() {
    let (mut sys_a, x) = mixed_system();
    let (mut sys_b, _) = mixed_system();

    // A: warm with compute_fg then call compute_f — cache hit expected.
    let mut g_a = vec![0.0; x.len()];
    let _ = compute_fg(&x, &mut sys_a, &mut g_a);
    let f_a = compute_f(&x, &mut sys_a);

    // B: always fresh.
    let mut g_b = vec![0.0; x.len()];
    force_cache_miss(&mut sys_b);
    let _ = compute_fg(&x, &mut sys_b, &mut g_b);
    force_cache_miss(&mut sys_b);
    let f_b = compute_f(&x, &mut sys_b);

    assert_eq!(f_a.to_bits(), f_b.to_bits());
    assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
    assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
}

// ── packer's "unscaled re-evaluation" pattern ─────────────────────────────
//
// After `pgencan` converges, `packer.rs` swaps `radius := radius_ini` and calls
// `compute_f` at the same `x` to measure violations under the true (unscaled)
// atomic radii. The cache key intentionally does not include `radius`, so this
// pattern hits the cache — verify the result under a radius mutation between
// the scaled and unscaled calls is identical to a fresh full evaluation.

#[test]
fn radii_swap_between_fg_and_f_cached_matches_fresh() {
    let (mut sys_a, x) = mixed_system();
    let (mut sys_b, _) = mixed_system();

    // Scaled radii (typical during packing: discale=1.2).
    let scaled: Vec<F> = sys_a.radius_ini.iter().map(|r| r * 1.2).collect();
    let unscaled = sys_a.radius_ini.clone();

    // --- A: cached path ---
    sys_a.radius = scaled.clone();
    sys_a.sync_atom_props();
    let mut g_a = vec![0.0; x.len()];
    let _ = compute_fg(&x, &mut sys_a, &mut g_a);

    sys_a.radius = unscaled.clone();
    sys_a.sync_atom_props();
    let f_a = compute_f(&x, &mut sys_a);

    // --- B: always rebuild ---
    sys_b.radius = scaled.clone();
    sys_b.sync_atom_props();
    let mut g_b = vec![0.0; x.len()];
    force_cache_miss(&mut sys_b);
    let _ = compute_fg(&x, &mut sys_b, &mut g_b);

    sys_b.radius = unscaled.clone();
    sys_b.sync_atom_props();
    force_cache_miss(&mut sys_b);
    let f_b = compute_f(&x, &mut sys_b);

    assert_eq!(f_a.to_bits(), f_b.to_bits());
    assert_eq!(sys_a.fdist.to_bits(), sys_b.fdist.to_bits());
    assert_eq!(sys_a.frest.to_bits(), sys_b.frest.to_bits());
}

// ── move_flag path must stay on the slow path ─────────────────────────────
//
// When `move_flag` is true (during movebad), per-atom `fdist_atom` /
// `frest_atom` are accumulated inside the pair / constraint kernels. Running
// the cache path a second time would double-count; so cache must be bypassed
// whenever `move_flag` is set.

#[test]
fn move_flag_true_bypasses_cache() {
    let (mut sys, x) = mixed_system();

    // Warm the cache at normal (move_flag=false) state.
    let _ = compute_f(&x, &mut sys);
    assert!(sys.work.cached_geometry_valid);

    // Now turn on move_flag and reset per-atom trackers.
    sys.move_flag = true;
    sys.fdist_atom.iter_mut().for_each(|v| *v = 0.0);
    sys.frest_atom.iter_mut().for_each(|v| *v = 0.0);
    let _ = compute_f(&x, &mut sys);
    let fdist_move_a = sys.fdist_atom.clone();
    let frest_move_a = sys.frest_atom.clone();

    // Fresh context, same sequence, cache forced off each call.
    let (mut sys2, _) = mixed_system();
    force_cache_miss(&mut sys2);
    let _ = compute_f(&x, &mut sys2);
    sys2.move_flag = true;
    sys2.fdist_atom.iter_mut().for_each(|v| *v = 0.0);
    sys2.frest_atom.iter_mut().for_each(|v| *v = 0.0);
    force_cache_miss(&mut sys2);
    let _ = compute_f(&x, &mut sys2);

    assert_eq!(
        fdist_move_a.len(),
        sys2.fdist_atom.len(),
        "fdist_atom shape mismatch"
    );
    for (i, (&a, &b)) in fdist_move_a.iter().zip(sys2.fdist_atom.iter()).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "fdist_atom[{i}] mismatch under move_flag: {a} vs {b}"
        );
    }
    for (i, (&a, &b)) in frest_move_a.iter().zip(sys2.frest_atom.iter()).enumerate() {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "frest_atom[{i}] mismatch under move_flag: {a} vs {b}"
        );
    }
}
