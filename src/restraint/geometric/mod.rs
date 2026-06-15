//! Concrete geometric restraint types (Packmol kinds 2–15).
//!
//! The 14 structs `impl AtomRestraint` with their value/gradient bodies moved
//! verbatim from the original single-file `src/restraint.rs`. They are split
//! across two leaf modules purely to stay within the per-file LOC budget:
//! [`bounded`] holds the cube/box/sphere/ellipsoid pairs (kinds 2–9) and
//! [`surface`] holds the plane/cylinder/Gaussian families (kinds 10–15). The
//! re-exports below keep every external path — `crate::restraint::<Type>` and
//! the crate-root `molpack::<Type>` — identical to the single-file layout.

mod bounded;
mod surface;

pub use bounded::{
    InsideBoxRestraint, InsideCubeRestraint, InsideEllipsoidRestraint, InsideSphereRestraint,
    OutsideBoxRestraint, OutsideCubeRestraint, OutsideEllipsoidRestraint, OutsideSphereRestraint,
};
pub use surface::{
    AboveGaussianRestraint, AbovePlaneRestraint, BelowGaussianRestraint, BelowPlaneRestraint,
    InsideCylinderRestraint, OutsideCylinderRestraint,
};
