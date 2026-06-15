//! Surface restraints (Packmol kinds 10–15): plane, cylinder, Gaussian.
//!
//! Moved verbatim from the original single-file `src/restraint.rs`; each
//! `impl AtomRestraint` resolves the trait through the unchanged
//! `crate::restraint` path.

use crate::restraint::AtomRestraint;
use molrs::types::F;

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

impl AtomRestraint for AbovePlaneRestraint {
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

impl AtomRestraint for BelowPlaneRestraint {
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

/// Packmol kind 12 — quadratic penalty forcing atom inside finite cylinder.
#[derive(Debug, Clone, Copy)]
pub struct InsideCylinderRestraint {
    pub center: [F; 3],
    pub axis: [F; 3],
    pub radius: F,
    pub length: F,
}

impl InsideCylinderRestraint {
    pub fn new(center: [F; 3], axis: [F; 3], radius: F, length: F) -> Self {
        // The axis is normalized by `|axis|` in `f`/`fg`; a zero vector yields
        // NaN. The contract is a non-zero axis direction.
        debug_assert!(
            axis[0] != 0.0 || axis[1] != 0.0 || axis[2] != 0.0,
            "cylinder axis must be a non-zero direction"
        );
        Self {
            center,
            axis,
            radius,
            length,
        }
    }
}

impl AtomRestraint for InsideCylinderRestraint {
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
        // The axis is normalized by `|axis|` in `f`/`fg`; a zero vector yields
        // NaN. The contract is a non-zero axis direction.
        debug_assert!(
            axis[0] != 0.0 || axis[1] != 0.0 || axis[2] != 0.0,
            "cylinder axis must be a non-zero direction"
        );
        Self {
            center,
            axis,
            radius,
            length,
        }
    }
}

impl AtomRestraint for OutsideCylinderRestraint {
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

/// Packmol kind 14 — quadratic penalty above Gaussian surface.
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

impl AtomRestraint for AboveGaussianRestraint {
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

/// Packmol kind 15 — quadratic penalty below Gaussian surface.
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

impl AtomRestraint for BelowGaussianRestraint {
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
