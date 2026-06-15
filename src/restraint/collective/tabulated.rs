//! Tabulated (arbitrary-prior) members of the distribution-matching family: the
//! reaction coordinate ξ is driven to **any** target density supplied as a grid
//! `{(ξᵢ, ρᵢ)}` — so a profile from theory or experiment (e.g. the exact
//! Gouy–Chapman counter-ion profile) is packable directly. Gaussian and
//! exponential are the parametric special cases of this general form.
//!
//! The quantile function `q(p) = F⁻¹(p)` is precomputed once from the grid: the
//! cumulative integral `F(ξ) = ∫ρ / ∫ρ` (trapezoid, normalised to `[0, 1]`) is
//! built on the grid points and inverted by linear interpolation. The shared
//! Wasserstein [`engine`](super::engine) then matches it like any other target.

use molrs::types::F;

use super::Restraint;
use super::geometry::{plane_match_f, plane_match_fg, point_match_f, point_match_fg};

/// Precomputed inverse-CDF (quantile) sampler for a tabulated target density.
#[derive(Debug, Clone)]
struct Quantile {
    /// Ascending grid of reaction-coordinate values.
    xi: Vec<F>,
    /// Normalised cumulative integral at each grid point (`cdf[0]=0`, `cdf[last]=1`).
    cdf: Vec<F>,
}

impl Quantile {
    /// Build from a density grid. `xs` must be strictly ascending with ≥ 2
    /// points, `rho` non-negative with positive total mass.
    ///
    /// # Panics
    /// On a too-short, non-ascending, negative, or zero-mass grid.
    fn from_grid(xs: &[F], rho: &[F]) -> Self {
        assert!(xs.len() >= 2, "tabulated grid needs at least 2 points");
        assert_eq!(xs.len(), rho.len(), "tabulated xs/rho length mismatch");
        assert!(
            xs.windows(2).all(|w| w[1] > w[0]),
            "tabulated grid xs must be strictly ascending"
        );
        assert!(rho.iter().all(|&r| r >= 0.0), "tabulated rho must be ≥ 0");

        let n = xs.len();
        let mut cdf = vec![0.0 as F; n];
        for i in 1..n {
            let dx = xs[i] - xs[i - 1];
            cdf[i] = cdf[i - 1] + 0.5 * (rho[i] + rho[i - 1]) * dx;
        }
        let total = cdf[n - 1];
        assert!(
            total > 0.0,
            "tabulated density must have positive total mass"
        );
        for c in cdf.iter_mut() {
            *c /= total;
        }
        Self {
            xi: xs.to_vec(),
            cdf,
        }
    }

    /// `q(p) = F⁻¹(p)` for `p ∈ (0, 1)` by linear interpolation on the CDF grid.
    #[inline]
    fn quantile(&self, p: F) -> F {
        let n = self.cdf.len();
        if p <= self.cdf[0] {
            return self.xi[0];
        }
        if p >= self.cdf[n - 1] {
            return self.xi[n - 1];
        }
        // First index whose CDF reaches `p` (CDF is non-decreasing).
        let i = self.cdf.partition_point(|&c| c < p);
        let (c0, c1) = (self.cdf[i - 1], self.cdf[i]);
        let (x0, x1) = (self.xi[i - 1], self.xi[i]);
        if c1 > c0 {
            x0 + (p - c0) / (c1 - c0) * (x1 - x0)
        } else {
            x0
        }
    }
}

// ============================================================================
// TabulatedPlane — arbitrary target density of the signed plane distance
// ============================================================================

/// Drive a species' signed distance to a plane, `ξ = x·n̂ − offset`, to an
/// arbitrary target density supplied as a grid `{(ξᵢ, ρᵢ)}`.
#[derive(Debug, Clone)]
pub struct TabulatedPlane {
    normal: [F; 3],
    offset: F,
    strength: F,
    quant: Quantile,
}

impl TabulatedPlane {
    /// - `normal` — plane normal (normalised internally; need not be unit).
    /// - `offset` — plane offset so `ξ = x·n̂ − offset`.
    /// - `strength` — overall penalty multiplier.
    /// - `xs`, `rho` — target density grid (ξ strictly ascending, ρ ≥ 0).
    ///
    /// # Panics
    /// If the normal is the zero vector or the grid is invalid.
    pub fn new(normal: [F; 3], offset: F, strength: F, xs: &[F], rho: &[F]) -> Self {
        let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        assert!(norm > 0.0, "TabulatedPlane normal must be non-zero");
        Self {
            normal: [normal[0] / norm, normal[1] / norm, normal[2] / norm],
            offset,
            strength,
            quant: Quantile::from_grid(xs, rho),
        }
    }
}

impl Restraint for TabulatedPlane {
    fn f(&self, coords: &[[F; 3]], _scale: F, _scale2: F) -> F {
        plane_match_f(coords, &self.normal, self.offset, self.strength, |p| {
            self.quant.quantile(p)
        })
    }

    fn fg(&self, coords: &[[F; 3]], _scale: F, _scale2: F, grads: &mut [[F; 3]]) -> F {
        plane_match_fg(
            coords,
            &self.normal,
            self.offset,
            self.strength,
            |p| self.quant.quantile(p),
            grads,
        )
    }

    fn name(&self) -> &'static str {
        "TabulatedPlane"
    }
}

// ============================================================================
// TabulatedPoint — arbitrary target density of the distance to a point
// ============================================================================

/// Drive a species' distance to a centre, `ξ = ‖x − c‖`, to an arbitrary target
/// density supplied as a grid `{(ξᵢ, ρᵢ)}` (a prescribed radial profile).
#[derive(Debug, Clone)]
pub struct TabulatedPoint {
    center: [F; 3],
    strength: F,
    quant: Quantile,
}

impl TabulatedPoint {
    /// - `center` — the point `c` distances are measured from.
    /// - `strength` — overall penalty multiplier.
    /// - `xs`, `rho` — target radial density grid (ξ = r strictly ascending, ρ ≥ 0).
    ///
    /// # Panics
    /// If the grid is invalid.
    pub fn new(center: [F; 3], strength: F, xs: &[F], rho: &[F]) -> Self {
        Self {
            center,
            strength,
            quant: Quantile::from_grid(xs, rho),
        }
    }
}

impl Restraint for TabulatedPoint {
    fn f(&self, coords: &[[F; 3]], _scale: F, _scale2: F) -> F {
        point_match_f(coords, &self.center, self.strength, |p| {
            self.quant.quantile(p)
        })
    }

    fn fg(&self, coords: &[[F; 3]], _scale: F, _scale2: F, grads: &mut [[F; 3]]) -> F {
        point_match_fg(
            coords,
            &self.center,
            self.strength,
            |p| self.quant.quantile(p),
            grads,
        )
    }

    fn name(&self) -> &'static str {
        "TabulatedPoint"
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{assert_fd_grad, rng_uniform};
    use super::*;

    /// A fine grid sampling a Gaussian density; the tabulated quantile should
    /// reproduce the closed-form Gaussian quantile to grid accuracy.
    fn gaussian_grid(mu: F, sigma: F) -> (Vec<F>, Vec<F>) {
        let lo = mu - 6.0 * sigma;
        let hi = mu + 6.0 * sigma;
        let n = 1200;
        let xs: Vec<F> = (0..n)
            .map(|i| lo + (hi - lo) * i as F / (n - 1) as F)
            .collect();
        let rho: Vec<F> = xs
            .iter()
            .map(|&x| (-0.5 * ((x - mu) / sigma).powi(2)).exp())
            .collect();
        (xs, rho)
    }

    #[test]
    fn quantile_recovers_gaussian() {
        use super::super::engine::probit;
        let (xs, rho) = gaussian_grid(20.0, 5.0);
        let q = Quantile::from_grid(&xs, &rho);
        for &p in &[0.1, 0.25, 0.5, 0.75, 0.9] {
            let want = 20.0 + 5.0 * probit(p);
            assert!(
                (q.quantile(p) - want).abs() < 0.05,
                "tabulated q({p})={}, want {want}",
                q.quantile(p)
            );
        }
    }

    #[test]
    fn plane_gradient_matches_finite_difference() {
        let (xs, rho) = gaussian_grid(20.0, 5.0);
        let r = TabulatedPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, &xs, &rho);
        let mut seed = 0x1357_2468u64;
        let coords: Vec<[F; 3]> = (0..20)
            .map(|i| {
                [
                    rng_uniform(&mut seed, 0.0, 10.0),
                    rng_uniform(&mut seed, 0.0, 10.0),
                    2.0 * i as F + rng_uniform(&mut seed, 0.0, 0.4),
                ]
            })
            .collect();
        assert_fd_grad(&r, &coords);
    }

    #[test]
    fn point_gradient_matches_finite_difference() {
        // A decaying radial target (exponential-like grid).
        let n = 400;
        let xs: Vec<F> = (0..n).map(|i| 0.05 * i as F).collect();
        let rho: Vec<F> = xs.iter().map(|&x| (-x / 4.0).exp()).collect();
        let r = TabulatedPoint::new([0.0, 0.0, 0.0], 1000.0, &xs, &rho);
        let mut seed = 0x2468_1357u64;
        let coords: Vec<[F; 3]> = (0..20)
            .map(|i| {
                let rad = 2.0 + 1.5 * i as F + rng_uniform(&mut seed, 0.0, 0.3);
                let t = rng_uniform(&mut seed, 0.3, 1.2);
                [
                    rad * t.cos(),
                    rad * t.sin(),
                    rng_uniform(&mut seed, -3.0, 3.0),
                ]
            })
            .collect();
        assert_fd_grad(&r, &coords);
    }

    #[test]
    fn penalty_zero_on_target_quantiles() {
        let (xs, rho) = gaussian_grid(20.0, 5.0);
        let r = TabulatedPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, &xs, &rho);
        let n = 50;
        let q = Quantile::from_grid(&xs, &rho);
        let coords: Vec<[F; 3]> = (0..n)
            .map(|k| [0.0, 0.0, q.quantile((k as F + 0.5) / n as F)])
            .collect();
        assert!(r.f(&coords, 1.0, 1.0) < 1e-6);
    }
}
