//! Evaluation mode + output types for the objective/restraint kernels.
//!
//! These were previously wrapped by a zero-sized `Constraints` facade; the
//! dispatch now lives directly in [`PackContext::evaluate`](super::PackContext::evaluate).

use molrs::types::F;

/// Evaluation mode for the objective/constraint kernels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalMode {
    /// Function value only.
    FOnly,
    /// Gradient only.
    GradientOnly,
    /// Function + gradient.
    FAndGradient,
    /// Restmol mode (same compute path as F+G, semantically explicit for callers).
    RestMol,
}

/// Unified evaluation output.
#[derive(Debug, Clone, Copy, Default)]
pub struct EvalOutput {
    pub f_total: F,
    pub fdist_max: F,
    pub frest_max: F,
}
