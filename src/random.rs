//! Precision-stable random helpers for packing paths.

use molrs::types::F;
use rand::{Rng, RngCore};

/// Draw a uniform random number in `[0, 1)` from an f64 stream, then cast to `F`.
///
/// `F` is currently `f64`, so the cast is a no-op; drawing from a fixed f64
/// stream regardless keeps the RNG trajectory stable if `F` is ever narrowed,
/// isolating true numeric-precision effects from type-dependent random draws.
#[inline]
pub fn uniform01(rng: &mut impl Rng) -> F {
    rng.random::<f64>() as F
}

/// Same as [`uniform01`], but for trait-object RNGs used by hook runners.
#[inline]
pub fn uniform01_core(rng: &mut dyn RngCore) -> F {
    let unit = (rng.next_u64() as f64) / ((u64::MAX as f64) + 1.0);
    unit as F
}
