//! Precision-stable random helpers for packing paths.
//!
//! molpack keeps its own random streams on purpose (module-responsibility
//! ruling 10): GENCAN's initial placement and perturbations draw in Packmol's
//! order, and growth / lattice / torsion-MC draw from their own separate
//! streams, so each algorithm's trajectory — and the Packmol-parity and
//! growth goldens that pin it bit for bit — never moves when another
//! algorithm changes how many numbers it draws.

use molrs::op::F;
use rand::Rng;
use rand::RngExt;

/// Draw a uniform random number in `[0, 1)` from an f64 stream, then cast to `F`.
///
/// `F` is currently `f64`, so the cast is a no-op; drawing from a fixed f64
/// stream regardless keeps the RNG trajectory stable if `F` is ever narrowed,
/// isolating true numeric-precision effects from type-dependent random draws.
#[inline]
pub fn uniform01(rng: &mut impl Rng) -> F {
    rng.random::<f64>() as F
}

/// A unit draw for trait-object RNGs, used by the in-loop torsion optimizer.
/// It maps `next_u64` onto `[0, 1)` directly, which is a different stream from
/// [`uniform01`]; the two stay separate so neither the GENCAN / growth
/// trajectories nor the torsion-MC trajectories change (bit parity, see the
/// module docs).
#[inline]
pub fn uniform01_core(rng: &mut dyn Rng) -> F {
    let unit = (rng.next_u64() as f64) / ((u64::MAX as f64) + 1.0);
    unit as F
}
