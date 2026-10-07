//! Collective (group-level) restraints — penalties that see **every copy of a
//! species at once**.
//!
//! Where a per-atom [`AtomRestraint`](crate::AtomRestraint) sees one atom at a
//! time and contributes an independent external field `∑ᵢ U(xᵢ)`, a
//! [`Restraint`] sees *every* copy of a species at once and returns a
//! single penalty whose gradient is **coupled across the whole group**.
//!
//! # Two families
//!
//! **1. Distribution matching.** The coupling is what lets a species *follow* a
//! target spatial distribution: a per-atom field built from a target density is
//! minimised by collapsing every atom onto the density's mode, whereas a
//! distribution-distance penalty is minimised when the empirical distribution
//! *equals* the target.
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
//! **2. Pairwise separation.** A distribution target says where copies should
//! be; it does not say how close two of them may come. A separation penalty
//! states a **lower bound on the distance between two copies of the same
//! species** and is silent otherwise — a local exclusion, not a global profile.
//! [`SelfSeparation`] is the member: it keeps a species' molecules from
//! clustering with each other. Its sites are per-copy centroids, produced by the
//! shared [`com`] reduction.
//!
//! The two families measure different things and compose: both may be attached
//! to the same species at once.
//!
//! # Evaluation context
//!
//! `f` / `fg` receive a [`GroupCtx`] alongside the coordinates. It carries what a
//! group-level term cannot recover from a flat coordinate slice: the packer's
//! two annealing scales, **how many atoms make one copy** (so `coords` can be
//! cut into molecules), and the **minimum-image convention** in force (so a term
//! that measures distances agrees with the pair loop across a periodic
//! boundary).
//!
//! **Gradient convention** mirrors [`AtomRestraint`](crate::AtomRestraint):
//! `fg` accumulates `∂L/∂coords[i]` INTO `grads[i]` with `+=`. `coords` and
//! `grads` have equal length, one entry per atom in the group, in the packer's
//! own order: **copy-major, atom-minor**.

use molrs::core::{Mic, SimBox};
use molrs::op::F;

// ============================================================================
// Evaluation context
// ============================================================================

/// What the packer knows at evaluation time and a group-level term cannot
/// recover from `coords` alone.
///
/// Captured once per objective evaluation and passed by value (it is `Copy` and
/// four words wide). Deliberately *not* stored on the restraint: the cell is
/// resolved after the targets are lowered, so a restraint that cached a box at
/// construction time could cache the wrong one.
#[derive(Debug, Clone, Copy)]
pub struct GroupCtx<'a> {
    /// Linear-penalty annealing scale (Packmol's two-scale contract).
    pub scale: F,
    /// Quadratic-penalty annealing scale.
    pub scale2: F,
    /// Atoms in one copy of this species. `coords` is copy-major, so
    /// `coords.chunks(natoms_per_copy)` yields one molecule at a time.
    pub natoms_per_copy: usize,
    /// The packing cell. A term that needs to *partition* space — a neighbour
    /// search over the group — builds its own [`CellGrid`] from this, the same
    /// primitive and the same box the pair loop uses.
    ///
    /// [`CellGrid`]: molrs::core::CellGrid
    pub cell: &'a SimBox,
    /// Minimum-image convention in force, derived from [`cell`](Self::cell).
    /// [`Mic::Free`] when no axis wraps.
    ///
    /// Redundant with `cell` on purpose: it is computed once per objective
    /// evaluation and shared, because deriving it costs two array builds on a
    /// triclinic box and every restraint in the group would otherwise repeat
    /// them. The cell stays the authority; this is a cached representation.
    pub mic: Mic,
}

// ============================================================================
// Trait
// ============================================================================

/// Group-level penalty over all copies of one species.
///
/// Unlike [`AtomRestraint`](crate::AtomRestraint), which is evaluated once
/// per atom with only that atom's coordinate, a `Restraint` is
/// evaluated once per group with the coordinates of *all* copies. Its gradient
/// may therefore couple every particle to every other — exactly what a
/// distribution-matching or separation penalty needs.
pub trait Restraint: Send + Sync + std::fmt::Debug {
    /// Penalty value for the group's current configuration.
    ///
    /// `coords[i]` is the Cartesian position of the `i`-th atom in the group.
    fn f(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>) -> F;

    /// Fused value + gradient. Accumulates `∂L/∂coords[i]` INTO `grads[i]`
    /// with `+=`; returns the same value `f` would. `grads.len() == coords.len()`.
    fn fg(&self, coords: &[[F; 3]], ctx: GroupCtx<'_>, grads: &mut [[F; 3]]) -> F;

    /// Does this restraint express a **bound** — a condition that is either met
    /// or not — rather than a target it only ever approaches?
    ///
    /// A bound's penalty is exactly zero once it holds, so it can be folded into
    /// the packer's restraint verdict `frest` and gate convergence: an
    /// unsatisfiable bound then shows up as a run that does not converge,
    /// instead of one that quietly ignores it.
    ///
    /// A distribution target is the opposite: the squared-Wasserstein penalty of
    /// a finite sample never reaches zero, so folding it into `frest` would mean
    /// a pack that reproduces its target profile to three decimals is reported
    /// as non-convergent. Hence the default is `false`, and only members of the
    /// separation family override it.
    fn is_bound(&self) -> bool {
        false
    }

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

mod com;
mod engine;
mod exponential;
mod gaussian;
mod geometry;
mod separation;
mod tabulated;

pub use exponential::{ExponentialPlane, ExponentialPoint};
pub use gaussian::{GaussianPlane, GaussianPoint};
pub use separation::SelfSeparation;
pub use tabulated::{TabulatedPlane, TabulatedPoint};

/// Shared test helpers for the concrete restraint types: a dependency-free RNG
/// and a finite-difference gradient check.
#[cfg(test)]
pub(super) mod test_fixtures {
    use super::{GroupCtx, Restraint};
    use molrs::core::{Mic, SimBox};
    use molrs::op::F;

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

    /// A free-boundary cube big enough to hold the test coordinates; the
    /// partitioning a restraint builds from it must not change any answer.
    pub(crate) fn free_box(side: F) -> SimBox {
        SimBox::cube(side, molrs::op::F3::zeros(3), [false; 3]).expect("test box")
    }

    /// Unit scales, free boundaries, `natoms_per_copy` atoms per copy.
    pub(crate) fn ctx_free(cell: &SimBox, natoms_per_copy: usize) -> GroupCtx<'_> {
        GroupCtx {
            scale: 1.0,
            scale2: 1.0,
            natoms_per_copy,
            cell,
            mic: Mic::Free,
        }
    }

    /// Central finite-difference check of the analytic gradient along every axis.
    /// Coordinates must be spread enough that a 1e-6 perturbation never crosses a
    /// rank swap (where the sorted-CDF objective is non-smooth).
    pub(crate) fn assert_fd_grad(r: &dyn Restraint, coords: &[[F; 3]]) {
        let cell = free_box(1_000.0);
        assert_fd_grad_in(r, coords, ctx_free(&cell, 1));
    }

    /// [`assert_fd_grad`] under a caller-supplied context.
    pub(crate) fn assert_fd_grad_in(r: &dyn Restraint, coords: &[[F; 3]], ctx: GroupCtx<'_>) {
        let mut analytic = vec![[0.0 as F; 3]; coords.len()];
        r.fg(coords, ctx, &mut analytic);
        let eps = 1e-6;
        for i in 0..coords.len() {
            for k in 0..3 {
                let mut plus = coords.to_vec();
                let mut minus = coords.to_vec();
                plus[i][k] += eps;
                minus[i][k] -= eps;
                let fd = (r.f(&plus, ctx) - r.f(&minus, ctx)) / (2.0 * eps);
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
