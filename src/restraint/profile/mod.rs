//! Profile-restraint family: a reaction coordinate ξ(x) plus the
//! distribution/penalty math and the composed [`ProfileRestraint`].
//!
//! Three leaf layers feed the composition: the [`Coordinate`] reaction
//! coordinate with its analytic gradient ([`coordinate`]); the closed-form
//! target-distribution penalties — Boltzmann inversion with shell-volume
//! Jacobian and density floor ([`distribution`], [`ProfilePenalty`]); and the
//! tabulated/spline profiles ([`spline`], [`TabulatedProfile`]). The composed
//! [`restraint`] layer pairs a coordinate with a [`ProfileTarget`] (analytic or
//! tabulated) and implements `Restraint` by the chain rule
//! `∇ₓU = (dU/dξ)·∇ξ`.

pub mod coordinate;
pub mod distribution;
pub mod restraint;
pub mod spline;

pub use coordinate::{Coordinate, CoordinateError, PbcWrap};
pub use distribution::{DensityFloor, Distribution, InputKind, ProfilePenalty, ShellJacobian};
pub use restraint::{ProfileRestraint, ProfileTarget};
pub use spline::TabulatedProfile;

#[cfg(test)]
mod tests;
