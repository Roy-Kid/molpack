//! Profile-restraint family: a reaction coordinate ξ(x) plus the
//! distribution/penalty math and (in a later sub-spec) the composed `Restraint`.
//!
//! Two leaf layers are exposed so far: the [`Coordinate`] reaction coordinate
//! with its analytic gradient ([`coordinate`]), and the closed-form
//! target-distribution penalties — Boltzmann inversion with shell-volume
//! Jacobian and density floor ([`distribution`], [`ProfilePenalty`]). The
//! tabulated/spline profiles (`-04`) and the composed `Restraint` (`-05`) are
//! later siblings.

pub mod coordinate;
pub mod distribution;
pub mod spline;

pub use coordinate::{Coordinate, CoordinateError, PbcWrap};
pub use distribution::{DensityFloor, Distribution, InputKind, ProfilePenalty, ShellJacobian};
pub use spline::TabulatedProfile;
