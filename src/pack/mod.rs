//! Rigid-body packing: initial placement, the bad-move heuristic, and
//! GENCAN descent.
//!
//! One family, peer of [`crate::grow`]. Placement calls the bad-move
//! heuristic, and both call [`restmol`]; the heuristic does not call
//! placement. Shared grid installation lives in [`crate::context::grid`],
//! so growth does not depend on this family.
//! Euler rotations and the optional in-loop optimizer stay outside: both
//! serve more than this driver.

pub(crate) mod gencan;
pub(crate) mod initial;
pub(crate) mod movebad;
pub(crate) mod restmol;

pub use gencan::gencan_pack::GencanPack;
