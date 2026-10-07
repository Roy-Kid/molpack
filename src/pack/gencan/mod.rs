//! GENCAN optimizer — faithful Rust port of `gencan.f` and `pgencan.f90`.
//!
//! Reference: Birgin & Martinez, Comp.Opt.Appl. 23:101-125, 2002.

use molrs::op::F;
mod cg;
pub(super) mod gencan_pack;
mod linesearch;
mod phases;
mod search;
mod solver;
mod spg;

use linesearch::TnLsScratch;
pub use search::gencan;

/// Stage name shared by [`solver::GencanStage`] and the phase step report.
pub(crate) const STAGE_NAME: &str = "gencan";

// ── Precision-aware floors shared by the GENCAN phases ─────────────────────
//
// Packmol calibrates its thresholds for double precision, which is the
// active precision here (`F = f64`). The `.max(eps)`-style floors are a
// defensive lower bound; under f64 they are no-ops, but they keep each
// threshold meaningful if `F` is ever narrowed.

/// The "effectively zero" level for an objective value or a squared
/// residual norm: `1e-10`, floored at `F::EPSILON`.
#[inline]
fn small_floor() -> F {
    (1.0e-10 as F).max(F::EPSILON)
}

/// The shortest norm still treated as non-zero: `√ε`.
#[inline]
fn near_zero_norm_floor() -> F {
    F::EPSILON.sqrt()
}

/// A divisor floor that only keeps a norm away from exact zero.
#[inline]
fn positive_norm_floor() -> F {
    F::MIN_POSITIVE
}

use crate::Objective;

/// Parameters for the GENCAN call (matches `easygencan` defaults from `pgencan.f90`).
pub struct GencanParams {
    pub epsgpsn: F,
    pub maxit: usize,
    pub maxfc: usize,
    pub delmin: F,
}

impl Default for GencanParams {
    fn default() -> Self {
        Self {
            epsgpsn: 1.0e-6,
            maxit: 20,
            maxfc: 200,     // 10 * maxit
            delmin: 1.0e-2, // Packmol easygencan default (gencan.f: delmin = 1.d-2)
        }
    }
}

/// Result of a GENCAN run.
pub struct GencanResult {
    pub f: F,
    /// Projected-gradient sup-norm and the iteration count. Read by the
    /// optimizer tests only; the phase loop decides convergence from
    /// `fdist` / `frest` instead.
    #[cfg_attr(not(test), allow(dead_code))]
    pub gpsupn: F,
    #[cfg_attr(not(test), allow(dead_code))]
    pub iter: usize,
    /// 0=converged(eucl), 1=converged(sup), 2=noFprogress, 3=noGprogress,
    /// 4=fSmall, 7=maxIter, 8=maxFeval, <0=error
    pub inform: i32,
}

/// Reusable work buffers for repeated `pgencan` calls.
pub struct GencanWorkspace {
    g: Vec<F>,
    ind: Vec<usize>,
    d: Vec<F>,
    s: Vec<F>,
    y: Vec<F>,
    cg_scratch: cg::CgScratch,
    spg_scratch: spg::SpgScratch,
    tnls_scratch: TnLsScratch,
}

impl GencanWorkspace {
    pub fn new() -> Self {
        Self {
            g: Vec::new(),
            ind: Vec::new(),
            d: Vec::new(),
            s: Vec::new(),
            y: Vec::new(),
            cg_scratch: cg::CgScratch::new(0),
            spg_scratch: spg::SpgScratch::new(0),
            tnls_scratch: TnLsScratch::new(0),
        }
    }

    fn ensure_len(&mut self, n: usize) {
        if self.g.len() != n {
            self.g.resize(n, 0.0);
        }
        if self.d.len() != n {
            self.d.resize(n, 0.0);
        }
        if self.s.len() != n {
            self.s.resize(n, 0.0);
        }
        if self.y.len() != n {
            self.y.resize(n, 0.0);
        }
        if self.ind.capacity() < n {
            self.ind.reserve(n - self.ind.capacity());
        }
    }
}

impl Default for GencanWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

/// Entry point — mirrors `pgencan.f90` → `easygencan` → `gencan`.
///
/// The variable layout is:
///   x[0..3N]   = COM positions (free molecules)
///   x[3N..6N]  = Euler angles (free molecules)
///
/// Bounds: COM variables are unbounded; Euler angles may be bounded by
/// `constrain_rotation` constraints (Packmol pgencan.f90).
pub fn pgencan(
    x: &mut [F],
    obj: &mut dyn Objective,
    params: &GencanParams,
    precision: F,
    workspace: &mut GencanWorkspace,
) -> GencanResult {
    let n = x.len();

    let mut l = vec![0.0 as F; n];
    let mut u = vec![0.0 as F; n];
    obj.bounds(&mut l, &mut u);

    gencan(x, &l, &u, obj, params, precision, workspace)
}

#[cfg(test)]
mod tests {
    //! Unit-level guards for the GENCAN optimizer (`src/pack/gencan/`).
    //!
    //! These tests pin the optimizer's behaviour on synthetic
    //! quadratics with known analytic minima — fast (< 50 ms each) and
    //! decoupled from the rest of the pack pipeline.
    //!
    //! The `Objective` trait is the only contract `gencan` depends on, so we
    //! provide a hand-rolled `Quadratic` impl rather than poking at
    //! `PackContext` internals. `fdist`/`frest` are reported as 0 so the
    //! Packmol-style early-exit check (`fdist < precision && frest < precision`)
    //! never fires when `precision = 0.0`; gencan then has to converge on its
    //! own gpsupn / maxit criterion.

    use crate::Objective;
    use crate::eval::{EvalMode, EvalOutput};
    use crate::pack::gencan::{GencanParams, GencanWorkspace, gencan, pgencan};
    use molrs::op::F;

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
                EvalMode::FAndGradient | EvalMode::GradientOnly => {
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

    // ── the seam markers the GENCAN stage declares ─────────────────────────────

    /// Owner-side half of acceptance ac-008: what `GencanStage` declares on the
    /// stage seam belongs here, not in `stage::tests` (which knows only fakes).
    ///
    /// `Placed::None` because the stage seeds its own placements with `initial()`
    /// — it needs nothing placed on entry; `Placed::All` because it returns with
    /// every free molecule placed. The two markers are what a pipeline checks
    /// before chaining anything after this stage.
    #[test]
    fn gencan_stage_requires_none_guarantees_all() {
        use crate::pack::gencan::solver::{GencanSettings, GencanStage};
        use crate::{Placed, Stage};

        // The declarations are construction-time constants: an empty system is
        // enough to read them, and using one keeps this test off the algorithm.
        let stage = GencanStage::new(GencanSettings::default(), Vec::new(), None, 0, 0);

        assert_eq!(
            stage.requires().placed,
            Placed::None,
            "GENCAN places the molecules itself, so it requires nothing placed"
        );
        assert_eq!(
            stage.guarantees().placed,
            Placed::All,
            "GENCAN returns with every free molecule placed"
        );
    }
}
