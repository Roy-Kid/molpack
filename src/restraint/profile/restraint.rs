//! Composed profile restraint: a reaction coordinate ξ(x) paired with a 1-D
//! distribution penalty U(ξ), implementing the `Restraint` trait by the chain
//! rule.
//!
//! A biased site assigned this restraint feels the penalty `U(ξ(x))` pulling it
//! toward the target profile. The gradient is the chain rule
//! `∇ₓU = (dU/dξ)·∇ξ(x)`, where ξ and ∇ξ come from the [`Coordinate`]
//! (sub-spec -02) and `U`, `dU/dξ` from the distribution side — either the
//! analytic [`ProfilePenalty`] (-03) or the [`TabulatedProfile`] (-04), unified
//! behind [`ProfileTarget`].
//!
//! **Weight choice.** The profile penalty `U` is a LINEAR-energy term (it has
//! energy units, not squared distance), so — matching Packmol's built-in linear
//! penalties (box / cube / plane) — it is multiplied through `scale`, NOT
//! `scale2` (which the quadratic overlap / region terms use). Riding `scale`
//! lets the bias follow the global radius soft-start so an early strong bias
//! never fights the overlap term; `ProfileRestraint` holds NO mutable ramp state.
//!
//! **Units.** `kt` is a FIXED energy-scale factor (not a simulated temperature)
//! sharing the overlap penalty's energy units; it sets the bias steepness.

use molrs::types::F;

use super::coordinate::{Coordinate, PbcWrap};
use super::distribution::ProfilePenalty;
use super::spline::TabulatedProfile;
use crate::restraint::Restraint;

/// Dispatch over the two distribution-side targets — the closed-form analytic
/// penalty and the tabulated spline — exposing a uniform `u` / `du_dxi` surface.
///
/// Named to avoid colliding with [`Distribution`](super::distribution::Distribution)
/// (the analytic-shape enum); both inner types already expose `u(ξ, kt)` and
/// `du_dxi(ξ, kt)`, so this enum only forwards.
#[derive(Debug, Clone)]
pub enum ProfileTarget {
    /// Closed-form Boltzmann-inversion penalty (Gaussian / erf / tanh /
    /// exponential, with shell Jacobian and density floor).
    Analytic(ProfilePenalty),
    /// Monotone cubic-spline profile of a tabulated target density.
    Tabulated(TabulatedProfile),
}

impl ProfileTarget {
    /// Penalty energy `U(ξ)` at reaction coordinate `xi` with energy scale `kt`.
    #[inline]
    pub fn u(&self, xi: F, kt: F) -> F {
        match self {
            ProfileTarget::Analytic(p) => p.u(xi, kt),
            ProfileTarget::Tabulated(t) => t.u(xi, kt),
        }
    }

    /// Force-generating derivative `dU/dξ` at `xi` with energy scale `kt`.
    #[inline]
    pub fn du_dxi(&self, xi: F, kt: F) -> F {
        match self {
            ProfileTarget::Analytic(p) => p.du_dxi(xi, kt),
            ProfileTarget::Tabulated(t) => t.du_dxi(xi, kt),
        }
    }
}

/// The single composed profile restraint: one [`Coordinate`] × one
/// [`ProfileTarget`], biasing a site toward a target 1-D distribution.
///
/// Immutable value type with no interior mutability, so it is `Send + Sync` and
/// safe to share read-only across parallel packing iterations. Derives `Debug`
/// (required by the [`Restraint`] trait bound).
///
/// See the [module docs](self) for the weight choice (linear `scale`) and units.
#[derive(Debug, Clone)]
pub struct ProfileRestraint {
    /// Reaction coordinate ξ(x) and its analytic gradient ∇ξ.
    coordinate: Coordinate,
    /// Distribution penalty U(ξ) and its derivative dU/dξ.
    target: ProfileTarget,
    /// Fixed energy-scale factor kT setting the bias steepness.
    kt: F,
    /// Minimum-image carrier for the coordinate's `x − reference` delta.
    pbc: PbcWrap,
}

impl ProfileRestraint {
    /// Compose a `coordinate` and a `target` into a restraint with energy scale
    /// `kt`, wrapping coordinate deltas with the given `pbc` carrier. Pass
    /// [`PbcWrap::none`] for a non-periodic system.
    pub fn new(coordinate: Coordinate, target: ProfileTarget, kt: F, pbc: PbcWrap) -> Self {
        Self {
            coordinate,
            target,
            kt,
            pbc,
        }
    }

    /// Ergonomic constructor for a non-periodic system: the coordinate's
    /// `x − reference` delta is the raw Cartesian difference (no wrap).
    pub fn non_periodic(coordinate: Coordinate, target: ProfileTarget, kt: F) -> Self {
        Self::new(coordinate, target, kt, PbcWrap::none())
    }
}

impl Restraint for ProfileRestraint {
    /// Value `scale · U(ξ(x))`. Multiplies through the LINEAR weight `scale`
    /// (profile U is a linear-energy term); `scale2` is unused. The coordinate's
    /// inner clamp and the distribution's density floor keep this finite by
    /// construction.
    fn f(&self, x: &[F; 3], scale: F, _scale2: F) -> F {
        let xi = self.coordinate.xi(x, &self.pbc);
        scale * self.target.u(xi, self.kt)
    }

    /// Fused value + gradient. Accumulates the chain-rule gradient
    /// `g[k] += scale · (dU/dξ) · ∇ξ[k]` INTO `g` with `+=` (never overwriting —
    /// many restraints may contribute to the same atom), then returns
    /// `self.f(...)` so the value is computed once per the trait contract.
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let xi = self.coordinate.xi(x, &self.pbc);
        let grad_xi = self.coordinate.grad_xi(x, &self.pbc);
        let du_dxi = self.target.du_dxi(xi, self.kt);
        for k in 0..3 {
            g[k] += scale * du_dxi * grad_xi[k];
        }
        self.f(x, scale, scale2)
    }

    /// Declare the periodic box when the coordinate is anchored in one (any
    /// active periodic axis in `pbc`), so the pair-kernel minimum-image wrap
    /// stays consistent. Mirrors [`InsideBoxRestraint`](crate::InsideBoxRestraint)'s
    /// return shape: `Some((min, max, periodic))` with `min = 0`, `max =
    /// pbc_length`. Defers to `None` for a non-periodic carrier.
    fn periodic_box(&self) -> Option<([F; 3], [F; 3], [bool; 3])> {
        if self.pbc.periodic.iter().any(|&p| p) {
            Some(([0.0; 3], self.pbc.pbc_length, self.pbc.periodic))
        } else {
            None
        }
    }
}
