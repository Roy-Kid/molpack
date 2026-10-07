//! Concrete geometric restraint types (Packmol kinds 2–15).
//!
//! The 14 structs `impl AtomRestraint` with their value/gradient bodies moved
//! verbatim from the original single-file `src/restraint.rs`. They are split
//! across two leaf modules purely to stay within the per-file LOC budget:
//! [`bounded`] holds the cube/box/sphere/ellipsoid pairs (kinds 2–9) and
//! [`surface`] holds the plane/cylinder families (kinds 10–13).
//!
//! Crate-private: these are the penalty kernels the `.inp` grammar lowers
//! onto (`script::build`). The public geometric restraint is
//! `RegionRestraint` over a molrs region.
//!
//! # Why these are not molrs regions
//!
//! Kept as Packmol-faithful kernels on purpose (module-responsibility ruling
//! 10): each body is the exact arithmetic of Packmol's `comprest` / `gwalls`
//! kind, operation for operation, so an `.inp` script reproduces Packmol's
//! penalty values — and therefore its packed coordinates — bit for bit. A
//! molrs region's signed distance is a different (if equivalent) formula and
//! would move those bits. Geometry questions that do not feed the penalty
//! value — such as whether a plane repeats along a cell vector — are asked of
//! molrs ([`HalfSpace`]).

mod bounded;
mod surface;

pub use bounded::{
    InsideBoxRestraint, InsideCubeRestraint, InsideEllipsoidRestraint, InsideSphereRestraint,
    OutsideBoxRestraint, OutsideCubeRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint,
};
pub use surface::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideCylinderRestraint, OutsideCylinderRestraint,
};

use molrs::core::{HalfSpace, Region};
use molrs::op::F;

/// Does the plane with normal `normal` repeat along `shift`? Asked of molrs's
/// [`HalfSpace::repeats_along`]; a zero normal bounds nothing, so it repeats
/// along every shift.
pub(crate) fn plane_repeats_along(normal: [F; 3], shift: [F; 3]) -> bool {
    HalfSpace::new(normal, [0.0; 3]).map_or(true, |plane| plane.repeats_along(shift))
}

#[cfg(test)]
mod tests {
    mod gradient;
    mod penalty;
}
