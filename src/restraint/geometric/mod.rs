//! Concrete geometric restraint types (Packmol kinds 2–15).
//!
//! The 14 structs `impl AtomRestraint` with their value/gradient bodies moved
//! verbatim from the original single-file `src/restraint.rs`. They are split
//! across two leaf modules purely to stay within the per-file LOC budget:
//! [`bounded`] holds the cube/box/sphere/ellipsoid pairs (kinds 2–9) and
//! [`surface`] holds the plane/cylinder families (kinds 10–13).
//!
//! Crate-private: these are the penalty kernels the `.inp` grammar lowers
//! onto (`script::build`) and the Packmol-example harness (`cases`) builds
//! directly. The public geometric restraint is `RegionRestraint` over a molrs
//! region.

mod bounded;
mod surface;

pub use bounded::{
    InsideBoxRestraint, InsideCubeRestraint, InsideEllipsoidRestraint, InsideSphereRestraint,
    OutsideBoxRestraint, OutsideCubeRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint,
};
pub use surface::{
    AbovePlaneRestraint, BelowPlaneRestraint, InsideCylinderRestraint, OutsideCylinderRestraint,
};

use molrs::types::F;

/// Does the half-space with unit outward `normal` repeat along `shift`?
///
/// True exactly when the shift runs parallel to the plane, so translating by
/// it moves no point across the boundary. The same test molrs's `HalfSpace`
/// applies to itself — these crate-private `.inp` kernels predate the region
/// lift and answer for themselves.
pub(crate) fn plane_repeats_along(normal: [F; 3], shift: [F; 3]) -> bool {
    let n2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
    let s2 = shift[0] * shift[0] + shift[1] * shift[1] + shift[2] * shift[2];
    if n2 == 0.0 || s2 == 0.0 {
        return true;
    }
    let crossing = normal[0] * shift[0] + normal[1] * shift[1] + normal[2] * shift[2];
    crossing.abs() <= 1e-9 * (n2 * s2).sqrt()
}

#[cfg(test)]
mod tests {
    mod gradient;
    mod penalty;
}
