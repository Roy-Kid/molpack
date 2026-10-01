//! Evaluation mode and the numbers one objective call returns.
//!
//! A leaf: the objective implements the call, and the context stores the
//! maxima. This module reaches into neither.

use molrs::types::F;

/// Evaluation mode for the shared objective.
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
    /// Objective value. Zero when the mode did not ask for it.
    pub f_total: F,
    /// Largest intermolecular contact violation left on the context.
    pub fdist_max: F,
    /// Largest restraint violation left on the context.
    pub frest_max: F,
}
