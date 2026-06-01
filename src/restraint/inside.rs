//! Inside-region restraints (cube / box / sphere / ellipsoid).
//!
//! Split from the original monolithic `restraint.rs`. See [`super`] for the
//! `Restraint` trait, the gradient/two-scale contract, and a custom-restraint
//! example.

use super::Restraint;
use molrs::types::F;

// ============================================================================
// Kind 2: InsideCubeRestraint
// ============================================================================

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

impl Restraint for InsideCubeRestraint {
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

// ============================================================================
// Kind 3: InsideBoxRestraint
// ============================================================================

/// Packmol kind 3 — quadratic penalty forcing atom inside axis-aligned box.
///
/// `periodic[k] == true` marks axis `k` as periodic: the pair-kernel
/// minimum-image wrap uses this box's extent for axis `k`, and the cell
/// list on axis `k` wraps instead of clamping. The soft confinement
/// penalty (`f` / `fg`) is always active regardless of `periodic` — a
/// periodic axis still keeps atoms in the reference image via the
/// penalty; `periodic` only affects pair-distance and cell lookup.
#[derive(Debug, Clone, Copy)]
pub struct InsideBoxRestraint {
    pub min: [F; 3],
    pub max: [F; 3],
    pub periodic: [bool; 3],
}

impl InsideBoxRestraint {
    /// Construct a box restraint with explicit bounds and per-axis
    /// periodicity flags. Pass `[false; 3]` for a purely-confining box;
    /// `[true; 3]` for a fully-periodic box; mixed flags give slab
    /// geometries (e.g. `[true, true, false]` = XY-periodic slab).
    pub fn new(min: [F; 3], max: [F; 3], periodic: [bool; 3]) -> Self {
        Self { min, max, periodic }
    }

    /// Construct an axis-aligned cubic box from an origin and isotropic
    /// side length, with per-axis periodicity flags.
    pub fn cube_from_origin(origin: [F; 3], side: F, periodic: [bool; 3]) -> Self {
        Self {
            min: origin,
            max: [origin[0] + side, origin[1] + side, origin[2] + side],
            periodic,
        }
    }

    /// Construct an axis-aligned box from a molrs-core `SimBox`.
    ///
    /// The restraint bounds are `origin` → `origin + lengths`; caller
    /// supplies the `periodic` flags explicitly (the `SimBox` type does
    /// not currently carry per-axis periodicity). Off-axis cells
    /// (triclinic tilts) are not representable as an axis-aligned
    /// restraint — write a custom `impl Restraint` for those.
    pub fn from_simbox(simbox: &molrs::region::SimBox, periodic: [bool; 3]) -> Self {
        let origin = simbox.origin_view();
        let lengths = simbox.lengths();
        let o = [origin[0], origin[1], origin[2]];
        Self {
            min: o,
            max: [o[0] + lengths[0], o[1] + lengths[1], o[2] + lengths[2]],
            periodic,
        }
    }
}

impl Restraint for InsideBoxRestraint {
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

    fn periodic_box(&self) -> Option<([F; 3], [F; 3], [bool; 3])> {
        if self.periodic.iter().any(|&p| p) {
            Some((self.min, self.max, self.periodic))
        } else {
            None
        }
    }
}

// ============================================================================
// Kind 4: InsideSphereRestraint
// ============================================================================

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

impl Restraint for InsideSphereRestraint {
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

// ============================================================================
// Kind 5: InsideEllipsoidRestraint
// ============================================================================

/// Packmol kind 5 — quadratic penalty forcing atom inside ellipsoid.
#[derive(Debug, Clone, Copy)]
pub struct InsideEllipsoidRestraint {
    pub center: [F; 3],
    pub axes: [F; 3],
    pub exponent: F,
}

impl InsideEllipsoidRestraint {
    pub fn new(center: [F; 3], axes: [F; 3], exponent: F) -> Self {
        Self {
            center,
            axes,
            exponent,
        }
    }
}

impl Restraint for InsideEllipsoidRestraint {
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
