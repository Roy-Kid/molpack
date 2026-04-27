//! Unit-level guards for the GENCAN optimizer (`src/gencan/`).
//!
//! Until now the entire optimizer was only exercised through the full
//! `Molpack::pack` driver, so a regression in CG / SPG / projected-gradient
//! handling could only surface as a numerical drift on the gated Packmol
//! suite. These tests pin the optimizer's behaviour on synthetic
//! quadratics with known analytic minima — fast (< 50 ms each) and
//! decoupled from the rest of the pack pipeline.
//!
//! The `Objective` trait is the only contract `gencan` depends on, so we
//! provide a hand-rolled `Quadratic` impl rather than poking at
//! `PackContext` internals. `fdist`/`frest` are reported as 0 so the
//! Packmol-style early-exit check (`fdist < precision && frest < precision`)
//! never fires when `precision = 0.0`; gencan then has to converge on its
//! own gpsupn / maxit criterion.

use molpack::constraints::{EvalMode, EvalOutput};
use molpack::gencan::{GencanParams, GencanWorkspace, gencan, pgencan};
use molpack::objective::Objective;
use molrs::types::F;

/// f(x) = 0.5 · Σ (xᵢ − μᵢ)²; ∇f = (x − μ); minimum at x = μ, f = 0.
struct Quadratic {
    mu: Vec<F>,
    ncf: usize,
    ncg: usize,
}

impl Quadratic {
    fn new(mu: Vec<F>) -> Self {
        Self { mu, ncf: 0, ncg: 0 }
    }
}

impl Objective for Quadratic {
    fn evaluate(&mut self, x: &[F], mode: EvalMode, gradient: Option<&mut [F]>) -> EvalOutput {
        debug_assert_eq!(x.len(), self.mu.len());

        let mut f = 0.0;
        for (xi, mu_i) in x.iter().zip(self.mu.iter()) {
            let dx = xi - mu_i;
            f += 0.5 * dx * dx;
        }

        match mode {
            EvalMode::FOnly => {
                self.ncf += 1;
            }
            EvalMode::FAndGradient | EvalMode::GradientOnly | EvalMode::RestMol => {
                self.ncf += 1;
                self.ncg += 1;
                if let Some(g) = gradient {
                    for ((gi, xi), mu_i) in g.iter_mut().zip(x.iter()).zip(self.mu.iter()) {
                        *gi = xi - mu_i;
                    }
                }
            }
        }

        EvalOutput {
            f_total: f,
            fdist_max: 0.0,
            frest_max: 0.0,
        }
    }

    fn fdist(&self) -> F {
        0.0
    }

    fn frest(&self) -> F {
        0.0
    }

    fn ncf(&self) -> usize {
        self.ncf
    }

    fn ncg(&self) -> usize {
        self.ncg
    }

    fn reset_eval_counters(&mut self) {
        self.ncf = 0;
        self.ncg = 0;
    }
}

// ── unconstrained CG path ──────────────────────────────────────────────────

/// Drive `pgencan` on a 3-D positive-definite quadratic with no bounds.
/// The CG inner solver should hit the analytic minimum (μ) within the
/// gpsupn tolerance after a handful of outer iterations. If a future
/// refactor breaks `cg_solve`'s descent direction or the outer
/// projected-gradient bookkeeping, x will not converge to μ.
#[test]
fn pgencan_unconstrained_quadratic_converges_to_minimum() {
    let mu = vec![3.0, -2.0, 1.5];
    let mut obj = Quadratic::new(mu.clone());
    let mut x = vec![0.0; mu.len()];
    let params = GencanParams::default();
    let mut ws = GencanWorkspace::new();

    // precision = 0.0 → packmolprecision (`fdist < 0 && frest < 0`) is
    // never satisfied; gencan terminates only on gpsupn/maxit.
    let result = pgencan(&mut x, &mut obj, &params, 0.0, &mut ws);

    assert!(
        result.inform == 0 || result.inform == 1,
        "gencan did not converge: inform={}, gpsupn={}, iter={}",
        result.inform,
        result.gpsupn,
        result.iter
    );
    for i in 0..x.len() {
        let err = (x[i] - mu[i]).abs();
        assert!(
            err < 1e-4,
            "x[{i}] = {} expected ≈ {} (err={err}, gpsupn={})",
            x[i],
            mu[i],
            result.gpsupn
        );
    }
    assert!(result.f < 1e-8, "f = {} expected ≈ 0", result.f);
}

// ── bound-constrained SPG path ─────────────────────────────────────────────

/// Drive `gencan` directly with explicit bounds so the projected
/// gradient hits a wall: minimum of f(x) = 0.5(x − 5)² over x ∈ [0, 1]
/// is at x = 1 (clamped, projected gradient = 0 at the upper bound).
/// Exercises the SPG / projection branch of gencan that the
/// unconstrained test above never reaches.
#[test]
fn gencan_bound_constrained_quadratic_clamps_to_boundary() {
    let mu = vec![5.0];
    let mut obj = Quadratic::new(mu.clone());
    let mut x = vec![0.0];
    let l = vec![0.0];
    let u = vec![1.0];
    let params = GencanParams::default();
    let mut ws = GencanWorkspace::new();

    let result = gencan(&mut x, &l, &u, &mut obj, &params, 0.0, &mut ws);

    assert!(
        result.inform >= 0,
        "gencan errored: inform={}",
        result.inform
    );
    let err = (x[0] - 1.0).abs();
    assert!(
        err < 1e-6,
        "expected x clamped to upper bound 1.0, got {} (err={err}, gpsupn={})",
        x[0],
        result.gpsupn
    );
    let expected_f = 0.5 * (1.0 - 5.0_f64).powi(2);
    assert!(
        (result.f - expected_f).abs() < 1e-6,
        "expected f = {expected_f}, got {}",
        result.f
    );
}

/// Bilateral clamp: minimum at μ = −3, search bounded on [0, 4] so the
/// unique constrained minimum is at the lower bound x = 0. Mirrors the
/// `_clamps_to_boundary` test on the opposite face.
#[test]
fn gencan_clamps_to_lower_bound() {
    let mut obj = Quadratic::new(vec![-3.0]);
    let mut x = vec![2.0];
    let l = vec![0.0];
    let u = vec![4.0];
    let params = GencanParams::default();
    let mut ws = GencanWorkspace::new();

    let result = gencan(&mut x, &l, &u, &mut obj, &params, 0.0, &mut ws);

    assert!(result.inform >= 0);
    assert!(
        (x[0] - 0.0).abs() < 1e-6,
        "expected x clamped to 0.0, got {}",
        x[0]
    );
}
