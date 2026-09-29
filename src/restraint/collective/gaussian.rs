//! Gaussian-target members of the distribution-matching family: the reaction
//! coordinate ξ is driven to a Gaussian distribution `𝒩(μ, σ²)` via the shared
//! Wasserstein engine. The target quantile function is `q(p) = μ + σ·Φ⁻¹(p)` for
//! both geometries; only ξ and its gradient (supplied by [`geometry`](super::geometry))
//! differ.

use molrs::types::F;

use super::engine::probit;
use super::geometry::{plane_match_f, plane_match_fg, point_match_f, point_match_fg, unit};
use super::{GroupCtx, Restraint};

/// Gaussian target quantile `q(p) = μ + σ·Φ⁻¹(p)`.
#[inline]
fn gaussian_quantile(p: F, mu: F, sigma: F) -> F {
    mu + sigma * probit(p)
}

// ============================================================================
// GaussianPlane — Gaussian distribution of the signed distance to a plane
// ============================================================================

/// Drive a species' signed distance to a plane, `ξ = x·n̂ − offset`, to a
/// Gaussian distribution `𝒩(μ, σ²)` — i.e. pack the copies into a **slab** of
/// centre `μ` and width `σ` along the plane normal.
///
/// `strength` (`λ`) scales the squared-Wasserstein penalty; the [`GroupCtx`] is
/// accepted for trait symmetry but unused (the target is fixed, not annealed
/// with the radius schedule).
#[derive(Debug, Clone)]
pub struct GaussianPlane {
    normal: [F; 3],
    offset: F,
    strength: F,
    mu: F,
    sigma: F,
}

impl GaussianPlane {
    /// - `normal` — plane normal (normalised internally; need not be unit).
    /// - `offset` — plane offset so `ξ = x·n̂ − offset`.
    /// - `strength` — overall multiplier `λ`.
    /// - `mu`, `sigma` — target Gaussian mean / standard deviation (`σ > 0`).
    ///
    /// # Panics
    /// If the normal is the zero vector or `sigma <= 0`.
    pub fn new(normal: [F; 3], offset: F, strength: F, mu: F, sigma: F) -> Self {
        assert!(sigma > 0.0, "GaussianPlane sigma must be positive");
        Self {
            normal: unit(normal),
            offset,
            strength,
            mu,
            sigma,
        }
    }

    #[inline]
    fn quantile(&self) -> impl Fn(F) -> F + '_ {
        move |p| gaussian_quantile(p, self.mu, self.sigma)
    }
}

impl Restraint for GaussianPlane {
    fn f(&self, coords: &[[F; 3]], _ctx: GroupCtx<'_>) -> F {
        plane_match_f(
            coords,
            &self.normal,
            self.offset,
            self.strength,
            self.quantile(),
        )
    }

    fn fg(&self, coords: &[[F; 3]], _ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F {
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
        "GaussianPlane"
    }
}

// ============================================================================
// GaussianPoint — Gaussian distribution of the distance to a point
// ============================================================================

/// Drive a species' distance to a centre, `ξ = ‖x − c‖`, to a Gaussian
/// distribution `𝒩(μ, σ²)` — i.e. pack the copies into a **spherical shell** of
/// radius `μ` and thickness `σ` around `center`.
///
/// This matches the distribution of the radius itself (the copies' radii are
/// `𝒩(μ, σ²)`-distributed). It is not a volumetric number density: a shell at
/// radius `μ` is what you get. Use `μ ≳ 3σ` so the Gaussian does not extend to
/// negative (unphysical) radii.
#[derive(Debug, Clone)]
pub struct GaussianPoint {
    center: [F; 3],
    strength: F,
    mu: F,
    sigma: F,
}

impl GaussianPoint {
    /// - `center` — the point `c` the shell is centred on.
    /// - `strength` — overall multiplier `λ`.
    /// - `mu`, `sigma` — target shell radius / thickness (`σ > 0`).
    ///
    /// # Panics
    /// If `sigma <= 0`.
    pub fn new(center: [F; 3], strength: F, mu: F, sigma: F) -> Self {
        assert!(sigma > 0.0, "GaussianPoint sigma must be positive");
        Self {
            center,
            strength,
            mu,
            sigma,
        }
    }

    #[inline]
    fn quantile(&self) -> impl Fn(F) -> F + '_ {
        move |p| gaussian_quantile(p, self.mu, self.sigma)
    }
}

impl Restraint for GaussianPoint {
    fn f(&self, coords: &[[F; 3]], _ctx: GroupCtx<'_>) -> F {
        point_match_f(coords, &self.center, self.strength, self.quantile())
    }

    fn fg(&self, coords: &[[F; 3]], _ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F {
        point_match_fg(coords, &self.center, self.strength, self.quantile(), grads)
    }

    fn name(&self) -> &'static str {
        "GaussianPoint"
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{assert_fd_grad, ctx_free, free_box, rng_uniform};

    #[test]
    fn a_distribution_target_is_not_a_bound() {
        // A finite sample's squared-Wasserstein penalty never reaches zero, so
        // folding it into `frest` would make a well-matched pack report as
        // non-convergent. Only the separation family opts in.
        assert!(!GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1.0, 0.0, 3.0).is_bound());
        assert!(!GaussianPoint::new([0.0; 3], 1.0, 5.0, 1.0).is_bound());
    }
    use super::*;

    #[test]
    fn plane_gradient_matches_finite_difference() {
        let r = GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, 20.0, 5.0);
        let mut seed = 0x1234_5678u64;
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
        let r = GaussianPoint::new([0.0, 0.0, 0.0], 1000.0, 30.0, 4.0);
        let mut seed = 0x9E37_79B9u64;
        let coords: Vec<[F; 3]> = (0..20)
            .map(|i| {
                let rad = 12.0 + 2.0 * i as F + rng_uniform(&mut seed, 0.0, 0.4);
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
    fn plane_penalty_zero_on_target_quantiles() {
        let r = GaussianPlane::new([0.0, 0.0, 1.0], 0.0, 1000.0, 20.0, 5.0);
        let n = 50;
        let coords: Vec<[F; 3]> = (0..n)
            .map(|k| [0.0, 0.0, 20.0 + 5.0 * probit((k as F + 0.5) / n as F)])
            .collect();
        assert!(r.f(&coords, ctx_free(&free_box(1_000.0), 1)) < 1e-6);
    }

    #[test]
    fn point_penalty_lower_for_shell_than_clump() {
        let r = GaussianPoint::new([0.0, 0.0, 0.0], 1000.0, 30.0, 4.0);
        let n = 60;
        let on_shell: Vec<[F; 3]> = (0..n)
            .map(|k| [30.0 + 4.0 * probit((k as F + 0.5) / n as F), 0.0, 0.0])
            .collect();
        let clump: Vec<[F; 3]> = (0..n).map(|_| [30.0, 0.0, 0.0]).collect();
        assert!(
            r.f(&on_shell, ctx_free(&free_box(1_000.0), 1))
                < r.f(&clump, ctx_free(&free_box(1_000.0), 1))
        );
    }
}
