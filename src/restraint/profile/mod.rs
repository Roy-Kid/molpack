//! Profile-restraint family: a reaction coordinate ξ(x) plus (in later
//! sub-specs) the distribution/penalty math and the composed `Restraint`.
//!
//! This module currently exposes only the leaf geometry layer — the
//! [`Coordinate`] reaction coordinate with its analytic gradient. Later siblings
//! append the penalty (`-03`/`-04`) and the composed restraint (`-05`).

pub mod coordinate;

pub use coordinate::{Coordinate, CoordinateError, PbcWrap};
