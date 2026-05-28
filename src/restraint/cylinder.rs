//! Finite-cylinder restraints (inside / outside).
//!
//! Split from the original monolithic `restraint.rs`. See [`super`] for the
//! `Restraint` trait, the gradient/two-scale contract, and a custom-restraint
//! example.

use super::Restraint;
use molrs::types::F;

// ============================================================================
// Kind 12: InsideCylinderRestraint
// ============================================================================

/// Packmol kind 12 — quadratic penalty confining an atom **inside** a finite
/// cylinder.
///
/// Let `p = pos - center`, `û = axis/‖axis‖`, axial coordinate `w = p·û`, and
/// squared radial distance `d = ‖p − w·û‖²`. The penalty sums three
/// independent quadratic violations (each active only when the atom is outside
/// that bound), scaled by `scale2`:
/// `scale2 · [ max(−w,0)² + max(w−length,0)² + max(d−radius²,0)² ]`
/// — below the base cap (`w<0`), past the top cap (`w>length`), or outside the
/// lateral wall (`d>radius²`).
///
/// - `center` — a point on the cylinder axis (the base-cap centre).
/// - `axis` — axis direction (need not be unit length; normalised internally).
/// - `radius` — cylinder radius.
/// - `length` — cylinder length measured along `axis` from `center`.
#[derive(Debug, Clone, Copy)]
pub struct InsideCylinderRestraint {
    pub center: [F; 3],
    pub axis: [F; 3],
    pub radius: F,
    pub length: F,
}

impl InsideCylinderRestraint {
    pub fn new(center: [F; 3], axis: [F; 3], radius: F, length: F) -> Self {
        Self {
            center,
            axis,
            radius,
            length,
        }
    }
}

impl Restraint for InsideCylinderRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axis;
        let (a1, a2, a3) = (x - c[0], y - c[1], z - c[2]);
        let vnorm = (ax[0].powi(2) + ax[1].powi(2) + ax[2].powi(2)).sqrt();
        let (vv1, vv2, vv3) = (ax[0] / vnorm, ax[1] / vnorm, ax[2] / vnorm);
        let w = vv1 * a1 + vv2 * a2 + vv3 * a3;
        let d = (a1 - vv1 * w).powi(2) + (a2 - vv2 * w).powi(2) + (a3 - vv3 * w).powi(2);
        scale2
            * ((-w).max(0.0).powi(2)
                + (w - self.length).max(0.0).powi(2)
                + (d - self.radius.powi(2)).max(0.0).powi(2))
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axis;
        let (a1, a2, a3) = (x - c[0], y - c[1], z - c[2]);
        let vnorm = (ax[0].powi(2) + ax[1].powi(2) + ax[2].powi(2)).sqrt();
        let (vv1, vv2, vv3) = (ax[0] / vnorm, ax[1] / vnorm, ax[2] / vnorm);
        let w = vv1 * a1 + vv2 * a2 + vv3 * a3;
        let d = (a1 - vv1 * w).powi(2) + (a2 - vv2 * w).powi(2) + (a3 - vv3 * w).powi(2);
        let len = self.length;
        let r2 = self.radius.powi(2);
        let rg0 = scale2
            * (-2.0 * (-w).max(0.0) * vv1
                + 2.0 * (w - len).max(0.0) * vv1
                + 2.0
                    * (d - r2).max(0.0)
                    * (2.0 * (a1 - vv1 * w) * (1.0 - vv1.powi(2))
                        + 2.0 * (a2 - vv2 * w) * (-vv2 * vv1)
                        + 2.0 * (a3 - vv3 * w) * (-vv3 * vv1)));
        let rg1 = scale2
            * (-2.0 * (-w).max(0.0) * vv2
                + 2.0 * (w - len).max(0.0) * vv2
                + 2.0
                    * (d - r2).max(0.0)
                    * (2.0 * (a1 - vv1 * w) * (-vv1 * vv2)
                        + 2.0 * (a2 - vv2 * w) * (1.0 - vv2.powi(2))
                        + 2.0 * (a3 - vv3 * w) * (-vv3 * vv2)));
        let rg2 = scale2
            * (-2.0 * (-w).max(0.0) * vv3
                + 2.0 * (w - len).max(0.0) * vv3
                + 2.0
                    * (d - r2).max(0.0)
                    * (2.0 * (a1 - vv1 * w) * (-vv1 * vv3)
                        + 2.0 * (a2 - vv2 * w) * (-vv2 * vv3)
                        + 2.0 * (a3 - vv3 * w) * (1.0 - vv3.powi(2))));
        g[0] += rg0;
        g[1] += rg1;
        g[2] += rg2;
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 13: OutsideCylinderRestraint
// ============================================================================

/// Packmol kind 13 — quadratic penalty forcing atom outside finite cylinder.
#[derive(Debug, Clone, Copy)]
pub struct OutsideCylinderRestraint {
    pub center: [F; 3],
    pub axis: [F; 3],
    pub radius: F,
    pub length: F,
}

impl OutsideCylinderRestraint {
    pub fn new(center: [F; 3], axis: [F; 3], radius: F, length: F) -> Self {
        Self {
            center,
            axis,
            radius,
            length,
        }
    }
}

impl Restraint for OutsideCylinderRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axis;
        let (a1, a2, a3) = (x - c[0], y - c[1], z - c[2]);
        let vnorm = (ax[0].powi(2) + ax[1].powi(2) + ax[2].powi(2)).sqrt();
        let (vv1, vv2, vv3) = (ax[0] / vnorm, ax[1] / vnorm, ax[2] / vnorm);
        let w = vv1 * a1 + vv2 * a2 + vv3 * a3;
        let d = (a1 - vv1 * w).powi(2) + (a2 - vv2 * w).powi(2) + (a3 - vv3 * w).powi(2);
        scale2
            * ((-w).min(0.0).powi(2)
                * (w - self.length).min(0.0).powi(2)
                * (d - self.radius.powi(2)).min(0.0).powi(2))
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axis;
        let (a1, a2, a3) = (x - c[0], y - c[1], z - c[2]);
        let vnorm = (ax[0].powi(2) + ax[1].powi(2) + ax[2].powi(2)).sqrt();
        let (vv1, vv2, vv3) = (ax[0] / vnorm, ax[1] / vnorm, ax[2] / vnorm);
        let w = vv1 * a1 + vv2 * a2 + vv3 * a3;
        let d = (a1 - vv1 * w).powi(2) + (a2 - vv2 * w).powi(2) + (a3 - vv3 * w).powi(2);
        let len = self.length;
        let r2 = self.radius.powi(2);
        let fra = (-w).min(0.0).powi(2);
        let frb = (w - len).min(0.0).powi(2);
        let frc = (d - r2).min(0.0).powi(2);
        let frab = fra * frb;
        let frac = fra * frc;
        let frbc = frb * frc;
        let dfra0 = -2.0 * (-w).min(0.0) * vv1;
        let dfrb0 = 2.0 * (w - len).min(0.0) * vv1;
        let dfrc0 = 2.0
            * (d - r2).min(0.0)
            * (2.0 * (a1 - vv1 * w) * (1.0 - vv1.powi(2))
                + 2.0 * (a2 - vv2 * w) * (-vv2 * vv1)
                + 2.0 * (a3 - vv3 * w) * (-vv3 * vv1));
        let dfra1 = -2.0 * (-w).min(0.0) * vv2;
        let dfrb1 = 2.0 * (w - len).min(0.0) * vv2;
        let dfrc1 = 2.0
            * (d - r2).min(0.0)
            * (2.0 * (a1 - vv1 * w) * (-vv1 * vv2)
                + 2.0 * (a2 - vv2 * w) * (1.0 - vv2.powi(2))
                + 2.0 * (a3 - vv3 * w) * (-vv3 * vv2));
        let dfra2 = -2.0 * (-w).min(0.0) * vv3;
        let dfrb2 = 2.0 * (w - len).min(0.0) * vv3;
        let dfrc2 = 2.0
            * (d - r2).min(0.0)
            * (2.0 * (a1 - vv1 * w) * (-vv1 * vv3)
                + 2.0 * (a2 - vv2 * w) * (-vv2 * vv3)
                + 2.0 * (a3 - vv3 * w) * (1.0 - vv3.powi(2)));
        g[0] += scale2 * (dfra0 * frbc + dfrb0 * frac + dfrc0 * frab);
        g[1] += scale2 * (dfra1 * frbc + dfrb1 * frac + dfrc1 * frab);
        g[2] += scale2 * (dfra2 * frbc + dfrb2 * frac + dfrc2 * frab);
        self.f(pos, scale, scale2)
    }
}
