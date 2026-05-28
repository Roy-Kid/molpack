//! Outside-region restraints (cube / box / sphere / ellipsoid).
//!
//! Split from the original monolithic `restraint.rs`. See [`super`] for the
//! `Restraint` trait, the gradient/two-scale contract, and a custom-restraint
//! example.

use super::Restraint;
use molrs::types::F;

// ============================================================================
// Kind 6: OutsideCubeRestraint
// ============================================================================

/// Packmol kind 6 — linear penalty forcing atom outside axis-aligned cube.
#[derive(Debug, Clone, Copy)]
pub struct OutsideCubeRestraint {
    pub origin: [F; 3],
    pub side: F,
}

impl OutsideCubeRestraint {
    pub fn new(origin: [F; 3], side: F) -> Self {
        Self { origin, side }
    }
}

impl Restraint for OutsideCubeRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.origin[0], self.origin[1], self.origin[2]);
        let (xmax, ymax, zmax) = (xmin + self.side, ymin + self.side, zmin + self.side);
        if x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax {
            let xmed = (xmax - xmin) / 2.0;
            let ymed = (ymax - ymin) / 2.0;
            let zmed = (zmax - zmin) / 2.0;
            let a1 = if x <= xmin + xmed { x - xmin } else { xmax - x };
            let a2 = if y <= ymin + ymed { y - ymin } else { ymax - y };
            let a3 = if z <= zmin + zmed { z - zmin } else { zmax - z };
            scale * (a1 + a2 + a3)
        } else {
            0.0
        }
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.origin[0], self.origin[1], self.origin[2]);
        let (xmax, ymax, zmax) = (xmin + self.side, ymin + self.side, zmin + self.side);
        if x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax {
            let xmed = (xmax - xmin) / 2.0;
            let ymed = (ymax - ymin) / 2.0;
            let zmed = (zmax - zmin) / 2.0;
            let (a1, a4) = if x <= xmin + xmed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            let (a2, a5) = if y <= ymin + ymed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            let (a3, a6) = if z <= zmin + zmed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            g[0] += scale * (a1 + a4);
            g[1] += scale * (a2 + a5);
            g[2] += scale * (a3 + a6);
        }
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 7: OutsideBoxRestraint
// ============================================================================

/// Packmol kind 7 — linear penalty forcing atom outside axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub struct OutsideBoxRestraint {
    pub min: [F; 3],
    pub max: [F; 3],
}

impl OutsideBoxRestraint {
    pub fn new(min: [F; 3], max: [F; 3]) -> Self {
        Self { min, max }
    }
}

impl Restraint for OutsideBoxRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.min[0], self.min[1], self.min[2]);
        let (xmax, ymax, zmax) = (self.max[0], self.max[1], self.max[2]);
        if x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax {
            let xmed = (xmax - xmin) / 2.0;
            let ymed = (ymax - ymin) / 2.0;
            let zmed = (zmax - zmin) / 2.0;
            let a1 = if x <= xmin + xmed { x - xmin } else { xmax - x };
            let a2 = if y <= ymin + ymed { y - ymin } else { ymax - y };
            let a3 = if z <= zmin + zmed { z - zmin } else { zmax - z };
            scale * (a1 + a2 + a3)
        } else {
            0.0
        }
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.min[0], self.min[1], self.min[2]);
        let (xmax, ymax, zmax) = (self.max[0], self.max[1], self.max[2]);
        if x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax {
            let xmed = (xmax - xmin) / 2.0;
            let ymed = (ymax - ymin) / 2.0;
            let zmed = (zmax - zmin) / 2.0;
            let (a1, a4) = if x <= xmin + xmed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            let (a2, a5) = if y <= ymin + ymed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            let (a3, a6) = if z <= zmin + zmed {
                (1.0, 0.0)
            } else {
                (0.0, -1.0)
            };
            g[0] += scale * (a1 + a4);
            g[1] += scale * (a2 + a5);
            g[2] += scale * (a3 + a6);
        }
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 8: OutsideSphereRestraint
// ============================================================================

/// Packmol kind 8 — quadratic penalty forcing atom outside sphere.
#[derive(Debug, Clone, Copy)]
pub struct OutsideSphereRestraint {
    pub center: [F; 3],
    pub radius: F,
}

impl OutsideSphereRestraint {
    pub fn new(center: [F; 3], radius: F) -> Self {
        Self { center, radius }
    }
}

impl Restraint for OutsideSphereRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let w = (x - c[0]).powi(2) + (y - c[1]).powi(2) + (z - c[2]).powi(2) - self.radius.powi(2);
        let a1 = w.min(0.0);
        scale2 * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let d = (x - c[0]).powi(2) + (y - c[1]).powi(2) + (z - c[2]).powi(2) - self.radius.powi(2);
        if d < 0.0 {
            g[0] += 4.0 * scale2 * (x - c[0]) * d;
            g[1] += 4.0 * scale2 * (y - c[1]) * d;
            g[2] += 4.0 * scale2 * (z - c[2]) * d;
        }
        self.f(pos, scale, scale2)
    }
}

// ============================================================================
// Kind 9: OutsideEllipsoidRestraint
// ============================================================================

/// Packmol kind 9 — quadratic penalty forcing atom outside ellipsoid.
///
/// Both `f` and `fg` apply the `scale2` factor — `f` previously omitted
/// it (a transcription artefact of the original Fortran-port comment),
/// which left `f` and `fg` 100× out of phase at the default
/// `scale2 = 0.01` and made the optimizer's gradient look 100× flatter
/// than the function value reported. The `scale2` factor here matches
/// the documented "quadratic penalty group" convention (kinds 4 / 5 /
/// 8 / 9 / 12 / 13) and the sister kind-5 [`InsideEllipsoidRestraint`].
#[derive(Debug, Clone, Copy)]
pub struct OutsideEllipsoidRestraint {
    pub center: [F; 3],
    pub axes: [F; 3],
    pub exponent: F,
}

impl OutsideEllipsoidRestraint {
    pub fn new(center: [F; 3], axes: [F; 3], exponent: F) -> Self {
        Self {
            center,
            axes,
            exponent,
        }
    }
}

impl Restraint for OutsideEllipsoidRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axes;
        let a1 = (x - c[0]).powi(2) / ax[0].powi(2);
        let a2 = (y - c[1]).powi(2) / ax[1].powi(2);
        let a3 = (z - c[2]).powi(2) / ax[2].powi(2);
        let w = a1 + a2 + a3 - self.exponent.powi(2);
        let v = w.min(0.0);
        scale2 * v * v
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let a1 = x - c[0];
        let b1 = y - c[1];
        let c1 = z - c[2];
        let a2 = self.axes[0].powi(2);
        let b2 = self.axes[1].powi(2);
        let c2 = self.axes[2].powi(2);
        let d = a1.powi(2) / a2 + b1.powi(2) / b2 + c1.powi(2) / c2 - self.exponent.powi(2);
        if d < 0.0 {
            let ds = scale2 * d;
            g[0] += 4.0 * ds * a1 / a2;
            g[1] += 4.0 * ds * b1 / b2;
            g[2] += 4.0 * ds * c1 / c2;
        }
        self.f(pos, scale, scale2)
    }
}
