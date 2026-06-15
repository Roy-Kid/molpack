//! Collective (group-level) restraints — the **geometric distribution-matching**
//! family.
//!
//! Where a per-atom [`AtomRestraint`](crate::restraint::AtomRestraint) sees one atom at a
//! time and contributes an independent external field `∑ᵢ U(xᵢ)`, a
//! [`Restraint`] sees *every* copy of a species at once and returns a
//! single penalty whose gradient is **coupled across the whole group**. That
//! coupling is what lets a species *follow* a target spatial distribution: a
//! per-atom field built from a target density is minimised by collapsing every
//! atom onto the density's mode, whereas a distribution-distance penalty is
//! minimised when the empirical distribution *equals* the target.
//!
//! # Structure of this family
//!
//! Every member matches a target **distribution** of a scalar reaction
//! coordinate ξ defined by a **geometry**, via the squared 1-D Wasserstein
//! (sorted-CDF) metric ([`engine`]). The two axes are orthogonal:
//!
//! - **geometry** ([`geometry`]) — maps Cartesian coordinates to ξ and scatters
//!   `∂L/∂ξ` back onto them: `plane` (ξ = signed distance to a plane → a slab),
//!   `point` (ξ = distance to a centre → a spherical shell), …
//! - **distribution** — the target quantile function `q(p) = F⁻¹(p)`: Gaussian,
//!   exponential, …
//!
//! Concrete types are the cross product, named `<Distribution><Geometry>`
//! and implementing [`Restraint`] directly (no wrapper, no builder —
//! same direction-3 convention as the per-atom restraints):
//! [`GaussianPlane`], [`GaussianPoint`], … Adding a distribution is a new
//! quantile function; adding a geometry is a new ξ/scatter pair; a new concrete
//! type then composes the two through the shared [`engine`].
//!
//! **Gradient convention** mirrors [`AtomRestraint`](crate::restraint::AtomRestraint):
//! `fg` accumulates `∂L/∂coords[i]` INTO `grads[i]` with `+=`. `coords` and
//! `grads` have equal length (one entry per atom in the group, same order).

use molrs::types::F;

// ============================================================================
// Trait
// ============================================================================

/// Group-level penalty over all copies of one species.
///
/// Unlike [`AtomRestraint`](crate::restraint::AtomRestraint), which is evaluated once
/// per atom with only that atom's coordinate, a `Restraint` is
/// evaluated once per group with the coordinates of *all* copies. Its gradient
/// may therefore couple every particle to every other — exactly what a
/// distribution-matching penalty needs.
pub trait Restraint: Send + Sync + std::fmt::Debug {
    /// Penalty value for the group's current configuration.
    ///
    /// `coords[i]` is the Cartesian position of the `i`-th atom in the group.
    fn f(&self, coords: &[[F; 3]], scale: F, scale2: F) -> F;

    /// Fused value + gradient. Accumulates `∂L/∂coords[i]` INTO `grads[i]`
    /// with `+=`; returns the same value `f` would. `grads.len() == coords.len()`.
    fn fg(&self, coords: &[[F; 3]], scale: F, scale2: F, grads: &mut [[F; 3]]) -> F;

    /// If `false`, the scheduler serializes this restraint (Python-backed
    /// collective restraints MUST return `false`).
    fn is_parallel_safe(&self) -> bool {
        true
    }

    /// Human-readable identifier.
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }
}

mod engine;
mod exponential;
mod gaussian;
mod geometry;
mod tabulated;

pub use exponential::{ExponentialPlane, ExponentialPoint};
pub use gaussian::{GaussianPlane, GaussianPoint};
pub use tabulated::{TabulatedPlane, TabulatedPoint};

/// Shared test helpers for the concrete `<Distribution><Geometry>` types:
/// a dependency-free RNG and a finite-difference gradient check.
#[cfg(test)]
pub(super) mod testutil {
    use super::Restraint;
    use molrs::types::F;

    /// Deterministic xorshift64* uniform in `[lo, hi)` — no external dep.
    pub(crate) fn rng_uniform(seed: &mut u64, lo: F, hi: F) -> F {
        let mut x = *seed;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        *seed = x;
        let u = (x.wrapping_mul(0x2545F4914F6CDD1D) >> 11) as F / (1u64 << 53) as F;
        lo + u * (hi - lo)
    }

    /// Central finite-difference check of the analytic gradient along every axis.
    /// Coordinates must be spread enough that a 1e-6 perturbation never crosses a
    /// rank swap (where the sorted-CDF objective is non-smooth).
    pub(crate) fn assert_fd_grad(r: &dyn Restraint, coords: &[[F; 3]]) {
        let mut analytic = vec![[0.0 as F; 3]; coords.len()];
        r.fg(coords, 1.0, 1.0, &mut analytic);
        let eps = 1e-6;
        for i in 0..coords.len() {
            for k in 0..3 {
                let mut plus = coords.to_vec();
                let mut minus = coords.to_vec();
                plus[i][k] += eps;
                minus[i][k] -= eps;
                let fd = (r.f(&plus, 1.0, 1.0) - r.f(&minus, 1.0, 1.0)) / (2.0 * eps);
                assert!(
                    (fd - analytic[i][k]).abs() < 1e-4,
                    "{} atom {i} axis {k}: fd={fd}, analytic={}",
                    r.name(),
                    analytic[i][k]
                );
            }
        }
    }
}
