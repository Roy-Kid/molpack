//! Shared distribution-matching engine: the squared 1-D Wasserstein (sorted-CDF)
//! penalty, plus the inverse standard-normal CDF (probit) used by Gaussian
//! targets.
//!
//! The penalty is geometry- and distribution-agnostic: it takes the group's
//! reaction coordinate `ξ` (already computed by a [`geometry`](super::geometry))
//! and a target quantile function `q(p) = F⁻¹(p)`. Sorting `ξ` and matching the
//! rank-`k` value to `q((k−½)/N)` yields
//!
//! ```text
//!   L = (λ / 2N) · ∑ₖ (ξ₍ₖ₎ − q_k)²,    ∂L/∂ξ₍ₖ₎ = (λ/N)(ξ₍ₖ₎ − q_k).
//! ```
//!
//! The assignment is recomputed every call (a sort), so the per-particle force
//! is a *linear* restoring force toward that particle's current target quantile
//! — coupled across the whole group, with fixed point `empirical distribution =
//! target`.

use molrs::types::F;

/// Value of the Wasserstein penalty for `xi` against the target whose quantile
/// function is `quantile`. `O(N log N)` (the sort).
pub(super) fn wasserstein_value(xi: &[F], strength: F, quantile: impl Fn(F) -> F) -> F {
    let n = xi.len();
    if n == 0 {
        return 0.0;
    }
    let order = sorted_order(xi);
    let inv_n = 1.0 / n as F;
    let mut sumsq = 0.0;
    for (k, &i) in order.iter().enumerate() {
        let d = xi[i] - quantile((k as F + 0.5) * inv_n);
        sumsq += d * d;
    }
    strength * 0.5 * sumsq * inv_n
}

/// Value + per-coordinate gradient of the Wasserstein penalty. Writes
/// `∂L/∂ξᵢ = (λ/N)(ξᵢ − q_{rank(i)})` into `dxi` (length `N`) and returns `L`.
/// A caller's geometry then projects `dxi` onto the Cartesian coordinates.
pub(super) fn wasserstein_grad(
    xi: &[F],
    strength: F,
    quantile: impl Fn(F) -> F,
    dxi: &mut [F],
) -> F {
    let n = xi.len();
    if n == 0 {
        return 0.0;
    }
    let order = sorted_order(xi);
    let inv_n = 1.0 / n as F;
    let mut sumsq = 0.0;
    for (k, &i) in order.iter().enumerate() {
        let d = xi[i] - quantile((k as F + 0.5) * inv_n);
        sumsq += d * d;
        dxi[i] = strength * inv_n * d;
    }
    strength * 0.5 * sumsq * inv_n
}

/// Indices of `xi` in ascending order (the rank → atom map).
#[inline]
fn sorted_order(xi: &[F]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..xi.len()).collect();
    order.sort_by(|&a, &b| xi[a].total_cmp(&xi[b]));
    order
}

/// Inverse standard-normal CDF (probit) via Acklam's rational approximation
/// (relative error < 1.2e-9). `p` must lie in `(0, 1)`.
pub(super) fn probit(p: F) -> F {
    const A: [F; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.38357751867269e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [F; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [F; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [F; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const PLOW: F = 0.02425;
    let phigh = 1.0 - PLOW;
    if p < PLOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= phigh {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probit_matches_known_quantiles() {
        assert!((probit(0.5)).abs() < 1e-9);
        assert!((probit(0.975) - 1.959963985).abs() < 1e-6);
        assert!((probit(0.025) + 1.959963985).abs() < 1e-6);
    }

    #[test]
    fn wasserstein_zero_when_on_quantiles() {
        // ξ placed exactly on the target quantiles → penalty ≈ 0.
        let n = 40;
        let q = |p: F| 10.0 + 3.0 * probit(p);
        let xi: Vec<F> = (0..n).map(|k| q((k as F + 0.5) / n as F)).collect();
        assert!(wasserstein_value(&xi, 1.0, q) < 1e-9);
    }
}
