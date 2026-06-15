//! Integration tests for the force-field geometry relaxer (`ff` feature).
//!
//! Gated on `ff`; without the feature the whole file compiles to an empty test
//! binary. Run with: `cargo test --features ff --test relaxer_lbfgs`.
#![cfg(feature = "ff")]

use std::sync::Arc;

use molpack::{F, InsideBoxRestraint, LBFGSRelaxer, Molpack, NullHandler, Potential, Target};

// ── A self-contained chain force field ───────────────────────────────────────

/// Harmonic bonds over consecutive atoms of a linear chain: `E = Σ ½k(r − r0)²`.
/// No angle or non-bonded terms, so the minimum is every consecutive bond at
/// exactly `r0` — easy to assert against after packing.
struct ChainHarmonic {
    n: usize,
    k: F,
    r0: F,
}

impl Potential for ChainHarmonic {
    fn calc_energy_forces(&self, coords: &[F]) -> (F, Vec<F>) {
        let mut e = 0.0;
        let mut f = vec![0.0; coords.len()];
        for i in 0..self.n - 1 {
            let (a, b) = (3 * i, 3 * (i + 1));
            let d = [
                coords[b] - coords[a],
                coords[b + 1] - coords[a + 1],
                coords[b + 2] - coords[a + 2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            e += 0.5 * self.k * (r - self.r0) * (r - self.r0);
            if r > 1e-12 {
                let coeff = self.k * (r - self.r0) / r; // dE/dr / r
                for t in 0..3 {
                    let fi = coeff * d[t];
                    f[a + t] += fi; // force on atom i   = -grad_i
                    f[b + t] -= fi; // force on atom i+1
                }
            }
        }
        (e, f)
    }
}

/// Zigzag chain coords with tetrahedral angles, `bond_length` apart.
fn chain_coords(n: usize, bond_length: F) -> Vec<[F; 3]> {
    let theta = 109.5 * std::f64::consts::PI as F / 180.0;
    let dx = bond_length * (-theta.cos());
    let dz = bond_length * theta.sin();
    let mut coords = Vec::with_capacity(n);
    coords.push([0.0, 0.0, 0.0]);
    for i in 1..n {
        let prev = coords[i - 1];
        let sign: F = if i % 2 == 0 { 1.0 } else { -1.0 };
        coords.push([prev[0] + dx, 0.0, prev[2] + sign * dz]);
    }
    coords
}

/// A chain force field over `n` atoms with the standard test stiffness.
fn chain_pot(n: usize, r0: F) -> Arc<dyn Potential> {
    Arc::new(ChainHarmonic { n, k: 100.0, r0 })
}

fn consecutive_distances(pos: &[[F; 3]]) -> Vec<F> {
    pos.windows(2)
        .map(|w| {
            let d = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .collect()
}

// ── tests ────────────────────────────────────────────────────────────────────

#[test]
fn target_accepts_forcefield_relaxer() {
    let n = 6;
    let coords = chain_coords(n, 1.5);
    let radii = vec![1.0; n];
    let pot = chain_pot(n, 1.2);
    // count == 1 is required (all copies share reference coords).
    let _t = Target::from_coords(&coords, &radii, 1).with_relaxer(LBFGSRelaxer::new(pot));
}

#[test]
#[should_panic(expected = "relaxers require count == 1")]
fn forcefield_relaxer_requires_count_1() {
    let n = 6;
    let coords = chain_coords(n, 1.5);
    let radii = vec![1.0; n];
    let pot = chain_pot(n, 1.2);
    let _t = Target::from_coords(&coords, &radii, 2).with_relaxer(LBFGSRelaxer::new(pot));
}

/// End-to-end: pack one flexible chain (force-field relaxer, `r0 = 1.2`)
/// alongside obstacle points so the packer must run iterations. The relaxer
/// fires each loop; shrinking the chain from its initial 1.5 Å bonds toward
/// `r0` never worsens packing, so the gate accepts it and the final geometry's
/// consecutive bonds reach `r0`. Rigid placement preserves bond lengths, so the
/// relaxation reaching the packed output proves the relaxer ran through the
/// packer loop and its result propagated to the result.
#[test]
fn pack_drives_forcefield_relaxation_to_r0() {
    let n = 8;
    let init_bond = 1.5;
    let r0 = 1.2;
    let coords = chain_coords(n, init_bond);
    let radii = vec![1.0; n];
    let pot = chain_pot(n, r0);

    let bx = InsideBoxRestraint::new([0.0, 0.0, 0.0], [12.0, 12.0, 12.0], [false; 3]);

    let chain = Target::from_coords(&coords, &radii, 1)
        .with_name("chain")
        .with_restraint(bx)
        .with_relaxer(LBFGSRelaxer::new(pot));

    // Obstacle points create packing pressure → the packer must iterate, so the
    // chain's relaxer fires inside the loop (a single small molecule in a roomy
    // box would converge at initialization and never enter the loop).
    let points = Target::from_coords(&[[0.0, 0.0, 0.0]], &[1.0], 30)
        .with_name("points")
        .with_restraint(bx);

    let result = Molpack::new()
        .with_handler(NullHandler)
        .with_seed(7)
        .pack_with_report(&[chain, points], 30)
        .expect("pack should not fail");

    assert_eq!(result.natoms(), n + 30);

    // The chain is the first target → its 8 atoms come first in the output.
    let chain_pos: Vec<[F; 3]> = result.positions()[..n].to_vec();
    let dists = consecutive_distances(&chain_pos);
    for (i, &d) in dists.iter().enumerate() {
        assert!(
            (d - r0).abs() < 0.05,
            "bond {i} should relax to r0={r0}, got {d:.4} (started at {init_bond})"
        );
    }
}
