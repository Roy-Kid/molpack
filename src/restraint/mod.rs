//! Restraint trait and concrete soft-penalty types for molecular packing.
//!
//! Each `*Restraint` struct is a **concrete, independent type** holding its own
//! geometric parameters — no `Builtin*` wrapper, no tagged-union `{kind, params[9]}`
//! blob, no builder pattern. User extensions `impl Restraint` identically and sit
//! beside the 14 Packmol-originals in type space.
//!
//! Numerical equivalence to the Fortran `comprest.f90` (value) and `gwalls.f90`
//! (gradient) is preserved branch-for-branch; see `docs/packmol_parity.md`.
//!
//! **Gradient convention**: `Restraint::fg` accumulates INTO `g` with `+=`.
//! Do not overwrite; many restraints may contribute to the same atom.
//!
//! **Two-scale contract** (Packmol convention): linear penalties
//! (box / cube / plane, kinds 2/3/6/7/10/11) use `scale`; quadratic penalties
//! (sphere / ellipsoid / cylinder / gaussian, kinds 4/5/8/9/12/13/14/15) use
//! `scale2`. Each `impl Restraint` decides internally which to consume.
//!
//! Direction-3 rule (see spec §0 bullet 9): all molrs-pack extension points
//! follow `pub trait X` + N concrete pub structs that `impl X`; user-defined
//! structs `impl X` the same way. No `Builtin*` prefix, no wrapper, no builder.
//!
//! # Writing your own restraint
//!
//! Implement [`Restraint`] on any `Debug` type. `f` returns the penalty value;
//! `fg` returns the same value AND accumulates the gradient INTO `g` with `+=`
//! (never `=`). Consume `scale` for linear penalties or `scale2` for quadratic
//! ones, matching the two-scale contract above.
//!
//! ```
//! use molpack::restraint::Restraint;
//!
//! /// Pull atoms toward the plane z = z0 with a quadratic well.
//! #[derive(Debug)]
//! struct ZWell { z0: f64 }
//!
//! impl Restraint for ZWell {
//!     fn f(&self, x: &[f64; 3], _scale: f64, scale2: f64) -> f64 {
//!         let dz = x[2] - self.z0;
//!         scale2 * dz * dz
//!     }
//!     fn fg(&self, x: &[f64; 3], scale: f64, scale2: f64, g: &mut [f64; 3]) -> f64 {
//!         let dz = x[2] - self.z0;
//!         g[2] += scale2 * 2.0 * dz; // += , not =
//!         self.f(x, scale, scale2)
//!     }
//! }
//! ```

use molrs::types::F;

mod cylinder;
mod gaussian;
mod inside;
mod outside;
mod plane;

pub use cylinder::{InsideCylinderRestraint, OutsideCylinderRestraint};
pub use gaussian::{AboveGaussianRestraint, BelowGaussianRestraint};
pub use inside::{
    InsideBoxRestraint, InsideCubeRestraint, InsideEllipsoidRestraint, InsideSphereRestraint,
};
pub use outside::{
    OutsideBoxRestraint, OutsideCubeRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint,
};
pub use plane::{AbovePlaneRestraint, BelowPlaneRestraint};

// ============================================================================
// Trait
// ============================================================================

/// Soft-penalty restraint evaluated per atom during packing.
///
/// - `f` — value only (line-search interpolation)
/// - `fg` — fused value + gradient; gradient accumulates INTO `g` with `+=`
/// - `is_parallel_safe` — if `false`, scheduler serializes this restraint
///   (Python-backed restraints MUST return `false`)
/// - `name` — human-readable identifier (default: `std::any::type_name::<Self>()`)
/// - `periodic_box` — opt-in: a restraint may declare that it defines a
///   periodic axis-aligned box (`min`, `max`, `periodic[k]` per axis).
///   At most one periodic box may be declared across all restraints on a
///   packing run; the packer resolves multiple declarations by requiring
///   identical bounds. Default `None` (non-periodic); only
///   [`InsideBoxRestraint`] overrides it.
///
/// `Debug` is required on concrete impls so `Target` / `Molpack` remain
/// printable for diagnostics. All built-in restraints derive it; user types
/// should do the same.
pub trait Restraint: Send + Sync + std::fmt::Debug {
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F;
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F;
    fn is_parallel_safe(&self) -> bool {
        true
    }
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }
    /// Declare that this restraint defines a periodic box for the
    /// pair-kernel minimum-image wrap. Return `Some((min, max, periodic))`
    /// where `periodic[k] == true` marks axis `k` as wrapping. Default
    /// `None` means this restraint does not imply PBC.
    fn periodic_box(&self) -> Option<([F; 3], [F; 3], [bool; 3])> {
        None
    }
}

/// Blanket impl so `Box<dyn Restraint>` itself implements the trait.
impl Restraint for Box<dyn Restraint> {
    #[inline]
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F {
        (**self).f(x, scale, scale2)
    }
    #[inline]
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        (**self).fg(x, scale, scale2, g)
    }
    #[inline]
    fn is_parallel_safe(&self) -> bool {
        (**self).is_parallel_safe()
    }
    #[inline]
    fn name(&self) -> &'static str {
        (**self).name()
    }
    #[inline]
    fn periodic_box(&self) -> Option<([F; 3], [F; 3], [bool; 3])> {
        (**self).periodic_box()
    }
}
