//! Half-space (plane) restraints.
//!
//! Split from the original monolithic `restraint.rs`. See [`super`] for the
//! `Restraint` trait, the gradient/two-scale contract, and a custom-restraint
//! example.

use super::Restraint;
use molrs::types::F;

// ============================================================================
// Kind 10: AbovePlaneRestraint
// ============================================================================

/// Packmol kind 10 — quadratic penalty forcing atom above plane `n·x >= d`.
#[derive(Debug, Clone, Copy)]
pub struct AbovePlaneRestraint {
    pub normal: [F; 3],
    pub distance: F,
}

impl AbovePlaneRestraint {
    pub fn new(normal: [F; 3], distance: F) -> Self {
        Self { normal, distance }
    }
}

impl Restraint for AbovePlaneRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let n = self.normal;
        let w = n[0] * pos[0] + n[1] * pos[1] + n[2] * pos[2] - self.distance;
        let a1 = w.min(0.0);
        scale * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let n = self.normal;
        let d = n[0] * pos[0] + n[1] * pos[1] + n[2] * pos[2] - self.distance;
        if d < 0.0 {
            let ds = scale * d;
            g[0] += 2.0 * n[0] * ds;
            g[1] += 2.0 * n[1] * ds;
            g[2] += 2.0 * n[2] * ds;
        }
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 11: BelowPlaneRestraint
// ============================================================================

/// Packmol kind 11 — quadratic penalty forcing atom below plane `n·x <= d`.
#[derive(Debug, Clone, Copy)]
pub struct BelowPlaneRestraint {
    pub normal: [F; 3],
    pub distance: F,
}

impl BelowPlaneRestraint {
    pub fn new(normal: [F; 3], distance: F) -> Self {
        Self { normal, distance }
    }
}

impl Restraint for BelowPlaneRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let n = self.normal;
        let w = n[0] * pos[0] + n[1] * pos[1] + n[2] * pos[2] - self.distance;
        let a1 = w.max(0.0);
        scale * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let n = self.normal;
        let d = n[0] * pos[0] + n[1] * pos[1] + n[2] * pos[2] - self.distance;
        if d > 0.0 {
            let ds = scale * d;
            g[0] += 2.0 * n[0] * ds;
            g[1] += 2.0 * n[1] * ds;
            g[2] += 2.0 * n[2] * ds;
        }
        self.f(pos, scale, scale2)
    }
}
