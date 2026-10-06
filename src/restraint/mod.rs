//! `AtomRestraint` trait and the soft penalties of molecular packing.
//!
//! Every public item here is re-exported at the crate root (one path per
//! item); the contracts a restraint implementer must keep are on
//! [`AtomRestraint`] itself.
//!
//! The `.inp` grammar's `inside box` / `outside sphere` / `above plane` …
//! keywords lower onto crate-private kernels in `geometric/` whose value and
//! gradient reproduce the Fortran `comprest.f90` / `gwalls.f90` branch for
//! branch (see `docs/packmol_parity.md`). They are the script layer's
//! implementation, not a second public vocabulary for shapes.
//!
//! Direction-3 rule: all molpack extension points
//! follow `pub trait X` + N concrete pub structs that `impl X`; user-defined
//! structs `impl X` the same way. No `Builtin*` prefix, no wrapper, no builder.

use molrs::op::types::F;
use molrs::spatial::SimBox;

// ============================================================================
// Trait
// ============================================================================

/// Soft-penalty restraint evaluated per atom during packing.
///
/// Geometry is not described by this trait. A region — a sphere, a box, a
/// cell, a mesh-bounded solid, a union of spheres, or any `&` / `|` / `~`
/// composition of them — is a molrs [`Region`](molrs::spatial::region::Region),
/// and the one public geometric restraint is [`RegionRestraint`]: stay inside
/// that region. [`CellRestraint`] is the same lift over a primitive cell plus
/// the lattice declaration the packer needs. User extensions `impl
/// AtomRestraint` for penalties that are not "stay inside a region".
///
/// **Gradient convention**: `fg` accumulates INTO `g` with `+=`. Do not
/// overwrite; many restraints may contribute to the same atom.
///
/// **Two-scale contract** (Packmol convention): linear penalties
/// (box / cube / plane, kinds 2/3/6/7/10/11) use `scale`; quadratic penalties
/// (sphere / ellipsoid / cylinder / gaussian, kinds 4/5/8/9/12/13/14/15) use
/// `scale2`. [`RegionRestraint`] is distance-quadratic and therefore in the
/// first class: it consumes `scale`. Each implementation decides internally
/// which to consume.
///
/// - `f` — value only (line-search interpolation)
/// - `fg` — fused value + gradient; gradient accumulates INTO `g` with `+=`
/// - `is_parallel_safe` — if `false`, scheduler serializes this restraint
///   (Python-backed restraints MUST return `false`)
/// - `name` — human-readable identifier (default: `std::any::type_name::<Self>()`)
/// - `declared_cell` — opt-in: a restraint confining atoms to a primitive
///   cell also states what that cell is. Periodicity itself is declared on
///   the engine entry (`with_periodic_box` / `with_cell`), never inferred
///   from a restraint's shape.
///
/// `Debug` is required on concrete impls so `Target` and the engines remain
/// printable for diagnostics. All built-in restraints derive it; user types
/// should do the same.
pub trait AtomRestraint: Send + Sync + std::fmt::Debug {
    fn f(&self, x: &[F; 3], scale: F, scale2: F) -> F;
    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F;
    fn is_parallel_safe(&self) -> bool {
        true
    }
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }
    /// Does this restraint still hold along `shift` when the cell wraps —
    /// that is, does it still mean the same thing in every image?
    ///
    /// Restraints are evaluated at lab coordinates, never at wrapped ones, so
    /// a restraint used with periodic boundaries has to answer the same way in
    /// every image. Two shapes do:
    ///
    /// - one that **confines** atoms inside a single image along `shift` — a
    ///   box, a sphere, a cell — because the atoms it allows never leave that
    ///   image;
    /// - one that **repeats** along `shift` — a half-space whose plane runs
    ///   parallel to it — because every image looks the same.
    ///
    /// Confining is only one of the two: a restraint that holds nothing in
    /// place still holds *along* a shift it repeats under.
    ///
    /// A restraint that does neither is open along `shift`: which side of it
    /// an atom falls on depends on which image the atom drifted into, and the
    /// packer refuses the combination by name rather than evaluating it in
    /// whichever image the origin happens to sit in.
    ///
    /// The default is `true`. A restraint molpack cannot inspect — a
    /// user-supplied `f` / `fg`, from Rust or from Python — is taken at its
    /// word here exactly as it is for every other property.
    fn holds_along(&self, shift: [F; 3]) -> bool {
        let _ = shift;
        true
    }

    /// The packing lattice this restraint defines, if it defines one.
    ///
    /// A restraint confining atoms to a primitive cell also states what that
    /// cell is, so the lattice and the confinement come from one declaration.
    fn declared_cell(&self) -> Option<SimBox> {
        None
    }
}

/// Blanket impl so `Box<dyn AtomRestraint>` itself implements the trait.
impl AtomRestraint for Box<dyn AtomRestraint> {
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
    fn declared_cell(&self) -> Option<SimBox> {
        (**self).declared_cell()
    }

    fn holds_along(&self, shift: [F; 3]) -> bool {
        (**self).holds_along(shift)
    }
}

pub(crate) mod cell;
mod collective;
pub(crate) mod geometric;
mod region;

pub use cell::CellRestraint;
pub use region::RegionRestraint;

pub use collective::{
    ExponentialPlane, ExponentialPoint, GaussianPlane, GaussianPoint, GroupCtx, Restraint,
    SelfSeparation, TabulatedPlane, TabulatedPoint,
};
