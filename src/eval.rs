//! Evaluation mode and the numbers one objective call returns.
//!
//! A leaf: the objective implements the call, and the system stores the
//! maxima. This module reaches into neither.

use molrs::op::F;

/// Evaluation mode for the shared objective.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalMode {
    /// Function value only.
    FOnly,
    /// Gradient only.
    GradientOnly,
    /// Function + gradient.
    FAndGradient,
}

/// Unified evaluation output.
#[derive(Debug, Clone, Copy, Default)]
pub struct EvalOutput {
    /// Objective value. Zero when the mode did not ask for it.
    pub f_total: F,
    /// Largest intermolecular contact violation left on the system.
    pub fdist_max: F,
    /// Largest restraint violation left on the system.
    pub frest_max: F,
}
