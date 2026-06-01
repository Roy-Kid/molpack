//! Gaussian-surface restraints (above / below).
//!
//! Split from the original monolithic `restraint.rs`. See [`super`] for the
//! `Restraint` trait, the gradient/two-scale contract, and a custom-restraint
//! example.

use super::Restraint;
use molrs::types::F;

// ============================================================================
// Kind 14: AboveGaussianRestraint
// ============================================================================

/// Packmol kind 14 — quadratic penalty keeping an atom **above** a Gaussian
/// bump surface `s(x,y) = z0 + height·exp(−(x−cx)²/2sx² − (y−cy)²/2sy²)`.
///
/// With signed gap `w = s(x,y) − z` (positive when the atom is below the
/// surface), the penalty is `scale · max(w,0)²` — active only while the atom
/// is below the surface. The exponent is floored at −50 to avoid underflow.
///
/// - `cx`, `cy` — surface peak position in x / y.
/// - `sx`, `sy` — Gaussian widths in x / y.
/// - `z0` — baseline z of the surface; `height` — bump amplitude.
#[derive(Debug, Clone, Copy)]
pub struct AboveGaussianRestraint {
    pub cx: F,
    pub cy: F,
    pub sx: F,
    pub sy: F,
    pub z0: F,
    pub height: F,
}

impl AboveGaussianRestraint {
    pub fn new(cx: F, cy: F, sx: F, sy: F, z0: F, height: F) -> Self {
        Self {
            cx,
            cy,
            sx,
            sy,
            z0,
            height,
        }
    }
}

impl Restraint for AboveGaussianRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let e1 = -(x - self.cx).powi(2) / (2.0 * self.sx.powi(2));
        let e2 = -(y - self.cy).powi(2) / (2.0 * self.sy.powi(2));
        let w = if e1 + e2 <= -50.0 {
            -(z - self.z0)
        } else {
            self.height * (e1 + e2).exp() - (z - self.z0)
        };
        let a1 = w.max(0.0);
        scale * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let e1 = -(x - self.cx).powi(2) / (2.0 * self.sx.powi(2));
        let e2 = -(y - self.cy).powi(2) / (2.0 * self.sy.powi(2));
        let d_raw = if e1 + e2 <= -50.0 {
            -(z - self.z0)
        } else {
            self.height * (e1 + e2).exp() - (z - self.z0)
        };
        if d_raw > 0.0 {
            let d = scale * d_raw;
            g[0] += -2.0 * d * (x - self.cx) * (d + (z - self.z0)) / self.sx.powi(2);
            g[1] += -2.0 * d * (y - self.cy) * (d + (z - self.z0)) / self.sy.powi(2);
            g[2] += -2.0 * d;
        }
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 15: BelowGaussianRestraint
// ============================================================================

/// Packmol kind 15 — quadratic penalty keeping an atom **below** a Gaussian
/// bump surface (same surface as [`AboveGaussianRestraint`]).
///
/// With signed gap `w = s(x,y) − z`, the penalty is `scale · min(w,0)²` —
/// active only while the atom is above the surface (`w < 0`). Parameters match
/// [`AboveGaussianRestraint`].
#[derive(Debug, Clone, Copy)]
pub struct BelowGaussianRestraint {
    pub cx: F,
    pub cy: F,
    pub sx: F,
    pub sy: F,
    pub z0: F,
    pub height: F,
}

impl BelowGaussianRestraint {
    pub fn new(cx: F, cy: F, sx: F, sy: F, z0: F, height: F) -> Self {
        Self {
            cx,
            cy,
            sx,
            sy,
            z0,
            height,
        }
    }
}

impl Restraint for BelowGaussianRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let e1 = -(x - self.cx).powi(2) / (2.0 * self.sx.powi(2));
        let e2 = -(y - self.cy).powi(2) / (2.0 * self.sy.powi(2));
        let w = if e1 + e2 <= -50.0 {
            -(z - self.z0)
        } else {
            self.height * (e1 + e2).exp() - (z - self.z0)
        };
        let a1 = w.min(0.0);
        scale * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let e1 = -(x - self.cx).powi(2) / (2.0 * self.sx.powi(2));
        let e2 = -(y - self.cy).powi(2) / (2.0 * self.sy.powi(2));
        let d_raw = if e1 + e2 <= -50.0 {
            -(z - self.z0)
        } else {
            self.height * (e1 + e2).exp() - (z - self.z0)
        };
        if d_raw < 0.0 {
            let d = scale * d_raw;
            g[0] += -2.0 * d * (x - self.cx) * (d + (z - self.z0)) / self.sx.powi(2);
            g[1] += -2.0 * d * (y - self.cy) * (d + (z - self.z0)) / self.sy.powi(2);
            g[2] += -2.0 * d;
        }
        self.f(pos, scale, scale2)
    }
}
