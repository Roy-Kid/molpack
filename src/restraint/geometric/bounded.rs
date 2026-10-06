//! Bounded-region restraints (Packmol kinds 2–9): cube, box, sphere, ellipsoid.
//!
//! Moved verbatim from the original single-file `src/restraint.rs`; each
//! `impl AtomRestraint` resolves the trait through the unchanged
//! `crate::restraint` path.

use crate::AtomRestraint;
use molrs::op::types::F;

/// Packmol kind 2 — quadratic penalty forcing atom inside axis-aligned cube.
#[derive(Debug, Clone, Copy)]
pub struct InsideCubeRestraint {
    pub origin: [F; 3],
    pub side: F,
}

impl InsideCubeRestraint {
    pub fn new(origin: [F; 3], side: F) -> Self {
        Self { origin, side }
    }
}

impl AtomRestraint for InsideCubeRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.origin[0], self.origin[1], self.origin[2]);
        let (xmax, ymax, zmax) = (xmin + self.side, ymin + self.side, zmin + self.side);
        let a1 = (x - xmin).min(0.0);
        let a2 = (y - ymin).min(0.0);
        let a3 = (z - zmin).min(0.0);
        let mut f = scale * (a1 * a1 + a2 * a2 + a3 * a3);
        let a1 = (x - xmax).max(0.0);
        let a2 = (y - ymax).max(0.0);
        let a3 = (z - zmax).max(0.0);
        f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
        f
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let (xmin, ymin, zmin) = (self.origin[0], self.origin[1], self.origin[2]);
        let (xmax, ymax, zmax) = (xmin + self.side, ymin + self.side, zmin + self.side);
        let a1 = x - xmin;
        let a2 = y - ymin;
        let a3 = z - zmin;
        if a1 < 0.0 {
            g[0] += scale * 2.0 * a1;
        }
        if a2 < 0.0 {
            g[1] += scale * 2.0 * a2;
        }
        if a3 < 0.0 {
            g[2] += scale * 2.0 * a3;
        }
        let a1 = x - xmax;
        let a2 = y - ymax;
        let a3 = z - zmax;
        if a1 > 0.0 {
            g[0] += scale * 2.0 * a1;
        }
        if a2 > 0.0 {
            g[1] += scale * 2.0 * a2;
        }
        if a3 > 0.0 {
            g[2] += scale * 2.0 * a3;
        }
        self.f(pos, scale, scale2)
    }
}

/// Packmol kind 3 — quadratic penalty forcing atom inside axis-aligned box.
///
/// Confinement only: periodicity is the engine's declaration
/// (`with_periodic_box` / `with_cell`), never a property of the shape that
/// confines. The penalty keeps atoms in the reference image either way.
#[derive(Debug, Clone, Copy)]
pub struct InsideBoxRestraint {
    pub min: [F; 3],
    pub max: [F; 3],
}

impl InsideBoxRestraint {
    /// Construct a box restraint with explicit bounds.
    pub fn new(min: [F; 3], max: [F; 3]) -> Self {
        Self { min, max }
    }
}

impl AtomRestraint for InsideBoxRestraint {
    fn f(&self, pos: &[F; 3], scale: F, _scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let a1 = (x - self.min[0]).min(0.0);
        let a2 = (y - self.min[1]).min(0.0);
        let a3 = (z - self.min[2]).min(0.0);
        let mut f = scale * (a1 * a1 + a2 * a2 + a3 * a3);
        let a1 = (x - self.max[0]).max(0.0);
        let a2 = (y - self.max[1]).max(0.0);
        let a3 = (z - self.max[2]).max(0.0);
        f += scale * (a1 * a1 + a2 * a2 + a3 * a3);
        f
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let a1 = x - self.min[0];
        let a2 = y - self.min[1];
        let a3 = z - self.min[2];
        if a1 < 0.0 {
            g[0] += scale * 2.0 * a1;
        }
        if a2 < 0.0 {
            g[1] += scale * 2.0 * a2;
        }
        if a3 < 0.0 {
            g[2] += scale * 2.0 * a3;
        }
        let a1 = x - self.max[0];
        let a2 = y - self.max[1];
        let a3 = z - self.max[2];
        if a1 > 0.0 {
            g[0] += scale * 2.0 * a1;
        }
        if a2 > 0.0 {
            g[1] += scale * 2.0 * a2;
        }
        if a3 > 0.0 {
            g[2] += scale * 2.0 * a3;
        }
        self.f(pos, scale, scale2)
    }
}

/// Packmol kind 4 — quadratic penalty forcing atom inside sphere.
#[derive(Debug, Clone, Copy)]
pub struct InsideSphereRestraint {
    pub center: [F; 3],
    pub radius: F,
}

impl InsideSphereRestraint {
    pub fn new(center: [F; 3], radius: F) -> Self {
        Self { center, radius }
    }
}

impl AtomRestraint for InsideSphereRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let w = (x - c[0]).powi(2) + (y - c[1]).powi(2) + (z - c[2]).powi(2) - self.radius.powi(2);
        let a1 = w.max(0.0);
        scale2 * a1 * a1
    }

    fn fg(&self, pos: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let d = (x - c[0]).powi(2) + (y - c[1]).powi(2) + (z - c[2]).powi(2) - self.radius.powi(2);
        if d > 0.0 {
            g[0] += 4.0 * scale2 * (x - c[0]) * d;
            g[1] += 4.0 * scale2 * (y - c[1]) * d;
            g[2] += 4.0 * scale2 * (z - c[2]) * d;
        }
        self.f(pos, scale, scale2)
    }
}

/// Packmol kind 5 — quadratic penalty forcing atom inside ellipsoid.
#[derive(Debug, Clone, Copy)]
pub struct InsideEllipsoidRestraint {
    pub center: [F; 3],
    pub axes: [F; 3],
    pub exponent: F,
}

impl InsideEllipsoidRestraint {
    pub fn new(center: [F; 3], axes: [F; 3], exponent: F) -> Self {
        // Each axis halves the penalty via division by `axes[k]²`; a zero axis
        // produces NaN/inf at eval time. The contract is positive semi-axes.
        debug_assert!(
            axes.iter().all(|&a| a > 0.0),
            "ellipsoid semi-axes must be strictly positive, got {axes:?}"
        );
        Self {
            center,
            axes,
            exponent,
        }
    }
}

impl AtomRestraint for InsideEllipsoidRestraint {
    fn f(&self, pos: &[F; 3], _scale: F, scale2: F) -> F {
        let (x, y, z) = (pos[0], pos[1], pos[2]);
        let c = self.center;
        let ax = self.axes;
        let a1 = (x - c[0]).powi(2) / ax[0].powi(2);
        let a2 = (y - c[1]).powi(2) / ax[1].powi(2);
        let a3 = (z - c[2]).powi(2) / ax[2].powi(2);
        let w = a1 + a2 + a3 - self.exponent.powi(2);
        let v = w.max(0.0);
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
        if d > 0.0 {
            g[0] += scale2 * 4.0 * d * a1 / a2;
            g[1] += scale2 * 4.0 * d * b1 / b2;
            g[2] += scale2 * 4.0 * d * c1 / c2;
        }
        self.f(pos, scale, scale2)
    }
}

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

impl AtomRestraint for OutsideCubeRestraint {
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

impl AtomRestraint for OutsideBoxRestraint {
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

impl AtomRestraint for OutsideSphereRestraint {
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
        // Each axis halves the penalty via division by `axes[k]²`; a zero axis
        // produces NaN/inf at eval time. The contract is positive semi-axes.
        debug_assert!(
            axes.iter().all(|&a| a > 0.0),
            "ellipsoid semi-axes must be strictly positive, got {axes:?}"
        );
        Self {
            center,
            axes,
            exponent,
        }
    }
}

impl AtomRestraint for OutsideEllipsoidRestraint {
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
