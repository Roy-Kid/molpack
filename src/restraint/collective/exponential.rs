//! Exponential-target members of the distribution-matching family: the reaction
//! coordinate ξ is driven to an exponential distribution (density `∝ e^(−ξ/λ)`,
//! `ξ ≥ 0`) via the shared Wasserstein engine — a smooth decay away from a wall
//! or centre (the diffuse-layer / Gouy–Chapman shape). The target quantile
//! function is `q(p) = −λ·ln(1 − p)` for both geometries; only ξ and its gradient
//! (supplied by [`geometry`](super::geometry)) differ.
//!
//! Because the quantiles are non-negative, both members place the densest copies
//! at the plane/centre (`ξ = 0`) and thin them out with decay length `λ`.

use molrs::types::F;

use super::Restraint;
use super::geometry::{plane_match_f, plane_match_fg, point_match_f, point_match_fg};

/// Exponential target quantile `q(p) = −λ·ln(1 − p)` (`≥ 0` for `p ∈ (0, 1)`).
#[inline]
fn exponential_quantile(p: F, lambda: F) -> F {
    -lambda * (1.0 - p).ln()
}

/// Validate / normalise a non-zero direction; panics on a zero vector.
fn unit(normal: [F; 3]) -> [F; 3] {
    let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    assert!(norm > 0.0, "normal must be non-zero");
    [normal[0] / norm, normal[1] / norm, normal[2] / norm]
}

// ============================================================================
// ExponentialPlane — exponential decay of the signed distance to a plane
// ============================================================================

/// Drive a species' signed distance to a plane, `ξ = x·n̂ − offset`, to an
/// exponential distribution `∝ e^(−ξ/λ)` (`ξ ≥ 0`) — a layer densest at the plane
/// and decaying along `+n̂` with decay length `λ` (a diffuse layer).
///
/// `strength` (`λ_pen`) scales the squared-Wasserstein penalty; `scale`/`scale2`
/// are accepted for trait symmetry but unused.
#[derive(Debug, Clone)]
pub struct ExponentialPlane {
    normal: [F; 3],
    offset: F,
    strength: F,
    lambda: F,
}

impl ExponentialPlane {
    /// - `normal` — plane normal (normalised internally; need not be unit).
    /// - `offset` — plane offset so `ξ = x·n̂ − offset` (the wall is at `ξ = 0`).
    /// - `strength` — overall penalty multiplier.
    /// - `lambda` — exponential decay length (`λ > 0`).
    ///
    /// # Panics
    /// If the normal is the zero vector or `lambda <= 0`.
    pub fn new(normal: [F; 3], offset: F, strength: F, lambda: F) -> Self {
        assert!(lambda > 0.0, "ExponentialPlane lambda must be positive");
        Self {
            normal: unit(normal),
            offset,
            strength,
            lambda,
        }
    }

    #[inline]
    fn quantile(&self) -> impl Fn(F) -> F + '_ {
        move |p| exponential_quantile(p, self.lambda)
    }
}

impl Restraint for ExponentialPlane {
    fn f(&self, coords: &[[F; 3]], _scale: F, _scale2: F) -> F {
        plane_match_f(
            coords,
            &self.normal,
            self.offset,
            self.strength,
            self.quantile(),
        )
    }

    fn fg(&self, coords: &[[F; 3]], _scale: F, _scale2: F, grads: &mut [[F; 3]]) -> F {
        plane_match_fg(
            coords,
            &self.normal,
            self.offset,
            self.strength,
            self.quantile(),
            grads,
        )
    }

    fn name(&self) -> &'static str {
        "ExponentialPlane"
    }
}

// ============================================================================
// ExponentialPoint — exponential decay of the distance to a point
// ============================================================================

/// Drive a species' distance to a centre, `ξ = ‖x − c‖`, to an exponential
/// distribution `∝ e^(−ξ/λ)` — densest at the centre and decaying radially with
/// decay length `λ` (a radial atmosphere around `center`).
#[derive(Debug, Clone)]
pub struct ExponentialPoint {
    center: [F; 3],
    strength: F,
    lambda: F,
}

impl ExponentialPoint {
    /// - `center` — the point `c` the decay is measured from.
    /// - `strength` — overall penalty multiplier.
    /// - `lambda` — radial decay length (`λ > 0`).
    ///
    /// # Panics
    /// If `lambda <= 0`.
    pub fn new(center: [F; 3], strength: F, lambda: F) -> Self {
        assert!(lambda > 0.0, "ExponentialPoint lambda must be positive");
        Self {
            center,
            strength,
            lambda,
        }
    }

    #[inline]
    fn quantile(&self) -> impl Fn(F) -> F + '_ {
        move |p| exponential_quantile(p, self.lambda)
    }
}

impl Restraint for ExponentialPoint {
    fn f(&self, coords: &[[F; 3]], _scale: F, _scale2: F) -> F {
        point_match_f(coords, &self.center, self.strength, self.quantile())
    }

    fn fg(&self, coords: &[[F; 3]], _scale: F, _scale2: F, grads: &mut [[F; 3]]) -> F {
        point_match_fg(coords, &self.center, self.strength, self.quantile(), grads)
    }

    fn name(&self) -> &'static str {
        "ExponentialPoint"
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{assert_fd_grad, rng_uniform};
    use super::*;

    #[test]
    fn quantile_is_nonnegative_and_increasing() {
        let q = |p| exponential_quantile(p, 3.0);
        assert!((q(0.0)).abs() < 1e-12); // q(0) = 0
        assert!(q(0.25) > 0.0 && q(0.75) > q(0.25));
    }

    #[test]
    fn plane_gradient_matches_finite_difference() {
        let r = ExponentialPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, 4.0);
        let mut seed = 0xABCD_1234u64;
        // Increasing, well-separated ξ so no rank swap under a 1e-6 perturbation.
        let coords: Vec<[F; 3]> = (0..20)
            .map(|i| {
                [
                    rng_uniform(&mut seed, 0.0, 8.0),
                    rng_uniform(&mut seed, 0.0, 8.0),
                    1.5 * i as F + rng_uniform(&mut seed, 0.0, 0.3),
                ]
            })
            .collect();
        assert_fd_grad(&r, &coords);
    }

    #[test]
    fn point_gradient_matches_finite_difference() {
        let r = ExponentialPoint::new([0.0, 0.0, 0.0], 1000.0, 4.0);
        let mut seed = 0x5151_7777u64;
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
        let r = ExponentialPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, 4.0);
        let n = 50;
        let coords: Vec<[F; 3]> = (0..n)
            .map(|k| [0.0, 0.0, exponential_quantile((k as F + 0.5) / n as F, 4.0)])
            .collect();
        assert!(r.f(&coords, 1.0, 1.0) < 1e-6);
    }
}
