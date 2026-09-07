//! Geometric `Region` trait and composition combinators.
//!
//! Mesh volumes: [`StlRegion`].
//!
//! A `Region` is a **geometric predicate** with a signed-distance function:
//! - `contains(x) == true`  ⇔ x is inside the region
//! - `signed_distance(x) < 0` inside, `> 0` outside, `== 0` on the boundary
//!
//! Regions compose into boolean combinations via `And` / `Or` / `Not`
//! (pure type algebra — zero runtime cost beyond the component evaluations).
//! Any `Region` lifts to a soft-penalty `AtomRestraint` via [`RegionRestraint`]:
//! `penalty(x) = scale2 * max(0, signed_distance(x))²`.
//!
//! # Example
//! ```
//! use molpack::region::{InsideSphereRegion, Not, Region, RegionExt};
//! use molpack::RegionRestraint;
//!
//! // Spherical shell: inside outer sphere AND NOT inside inner sphere
//! let shell = InsideSphereRegion::new([0.0; 3], 10.0)
//!     .and(Not(InsideSphereRegion::new([0.0; 3], 5.0)));
//! assert!(shell.contains(&[7.0, 0.0, 0.0]));   // in the shell
//! assert!(!shell.contains(&[3.0, 0.0, 0.0]));  // inside inner sphere
//! assert!(!shell.contains(&[15.0, 0.0, 0.0])); // outside outer sphere
//!
//! // Lift to a AtomRestraint for packing
//! let restraint = RegionRestraint(shell);
//! ```
//!
//! Direction-3 rule (spec §0 bullet 9): `Region` is a **separate trait**
//! from `AtomRestraint`; composition operators (`.and()` etc.) live on `Region`,
//! not on `AtomRestraint`. Plugin vs built-in `Region` are type-equal via
//! user `impl Region`.

use molrs::spatial::simbox::SimBox;
use molrs::types::F;
use ndarray::array;

use crate::restraint::AtomRestraint;

mod bvh;
mod stl;
pub use stl::{StlError, StlRegion};

// ============================================================================
// Core trait
// ============================================================================

/// Axis-aligned bounding box (AABB), used as an optimization hint for
/// cell-list setup. Optional (default `None`).
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: [F; 3],
    pub max: [F; 3],
}

/// A lattice declaration: `(H, origin, pbc)`, with the lattice vectors as the
/// columns of `H`.
pub type CellDeclaration = ([[F; 3]; 3], [F; 3], [bool; 3]);

/// Geometric predicate with signed-distance function.
pub trait Region: Send + Sync + std::fmt::Debug {
    /// Membership test. Must be consistent with `signed_distance(x) <= 0`.
    fn contains(&self, x: &[F; 3]) -> bool;

    /// Signed distance to the region boundary.
    /// - Negative inside, positive outside, zero on the boundary.
    /// - Not required to be the *Euclidean* signed distance for arbitrary
    ///   regions; only the sign and gradient direction matter for packing.
    fn signed_distance(&self, x: &[F; 3]) -> F;

    /// Gradient of `signed_distance` at `x`.
    ///
    /// Default implementation is a 3-point central finite difference
    /// (ε = 1e-6). Concrete `Region` types should override with an
    /// analytic gradient for hot-path use.
    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let h: F = 1e-6;
        let mut g = [0.0; 3];
        for k in 0..3 {
            let mut xp = *x;
            xp[k] += h;
            let mut xm = *x;
            xm[k] -= h;
            g[k] = (self.signed_distance(&xp) - self.signed_distance(&xm)) / (2.0 * h);
        }
        g
    }

    /// The packing lattice this region defines, if it defines one.
    ///
    /// A region confined to a primitive cell also *is* the cell: returning it
    /// here lets the packer pick up the lattice from the same declaration that
    /// confines the molecules, instead of making the caller state it twice and
    /// keeping the two in sync. Returns `(H, origin, pbc)` with the lattice
    /// vectors as columns of `H`.
    fn declared_cell(&self) -> Option<CellDeclaration> {
        None
    }

    /// Axis-aligned bounding box of the region. Used as an initialization
    /// hint; default `None` is always safe.
    fn bounding_box(&self) -> Option<Aabb> {
        None
    }

    /// True when this region is (or composes) a closed triangle mesh.
    /// Default `false`.
    fn is_closed_mesh(&self) -> bool {
        false
    }
}

// ============================================================================
// Combinators
// ============================================================================

/// Intersection of two regions: inside iff BOTH are inside.
/// Signed distance uses `max` (outside dominates — chain-rule selects the
/// larger component's gradient).
#[derive(Debug, Clone, Copy)]
pub struct And<A, B>(pub A, pub B);

impl<A: Region, B: Region> Region for And<A, B> {
    fn contains(&self, x: &[F; 3]) -> bool {
        self.0.contains(x) && self.1.contains(x)
    }
    fn signed_distance(&self, x: &[F; 3]) -> F {
        self.0.signed_distance(x).max(self.1.signed_distance(x))
    }
    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        if self.0.signed_distance(x) >= self.1.signed_distance(x) {
            self.0.signed_distance_grad(x)
        } else {
            self.1.signed_distance_grad(x)
        }
    }
    fn is_closed_mesh(&self) -> bool {
        self.0.is_closed_mesh() || self.1.is_closed_mesh()
    }
}

/// Union of two regions: inside iff EITHER is inside.
/// Signed distance uses `min` (inside dominates).
#[derive(Debug, Clone, Copy)]
pub struct Or<A, B>(pub A, pub B);

impl<A: Region, B: Region> Region for Or<A, B> {
    fn contains(&self, x: &[F; 3]) -> bool {
        self.0.contains(x) || self.1.contains(x)
    }
    fn signed_distance(&self, x: &[F; 3]) -> F {
        self.0.signed_distance(x).min(self.1.signed_distance(x))
    }
    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        if self.0.signed_distance(x) <= self.1.signed_distance(x) {
            self.0.signed_distance_grad(x)
        } else {
            self.1.signed_distance_grad(x)
        }
    }
    fn is_closed_mesh(&self) -> bool {
        self.0.is_closed_mesh() || self.1.is_closed_mesh()
    }
}

/// Complement: inside iff the inner region is NOT inside.
/// Signed distance is negated.
#[derive(Debug, Clone, Copy)]
pub struct Not<A>(pub A);

impl<A: Region> Region for Not<A> {
    fn contains(&self, x: &[F; 3]) -> bool {
        !self.0.contains(x)
    }
    fn signed_distance(&self, x: &[F; 3]) -> F {
        -self.0.signed_distance(x)
    }
    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let g = self.0.signed_distance_grad(x);
        [-g[0], -g[1], -g[2]]
    }
    fn is_closed_mesh(&self) -> bool {
        self.0.is_closed_mesh()
    }
}

// ============================================================================
// Extension trait for `.and() / .or() / .not()` chaining
// ============================================================================

/// Ergonomic method-chaining for `Region`.
///
/// ```
/// use molpack::region::{InsideSphereRegion, RegionExt};
/// let outer = InsideSphereRegion::new([0.0; 3], 10.0);
/// let inner = InsideSphereRegion::new([0.0; 3], 5.0);
/// let shell = outer.and(inner.not());
/// # let _ = shell;
/// ```
pub trait RegionExt: Region + Sized {
    fn and<B: Region>(self, other: B) -> And<Self, B> {
        And(self, other)
    }
    fn or<B: Region>(self, other: B) -> Or<Self, B> {
        Or(self, other)
    }
    fn not(self) -> Not<Self> {
        Not(self)
    }

    /// Lift this region into a soft-penalty [`AtomRestraint`]. Equivalent
    /// to wrapping in [`RegionRestraint`] manually.
    fn into_restraint(self) -> RegionRestraint<Self> {
        RegionRestraint(self)
    }
}

impl<R: Region + Sized> RegionExt for R {}

// ============================================================================
// RegionRestraint — lift any Region to a quadratic-exterior-penalty AtomRestraint
// ============================================================================

/// Wraps a `Region` as a soft-penalty `AtomRestraint` with quadratic
/// exterior penalty:
///
/// ```text
/// penalty(x) = scale2 * max(0, signed_distance(x))²
/// ```
///
/// Gradient uses the analytic chain rule
/// `2 * scale2 * max(0, d) * ∂d/∂x`, where `∂d/∂x` comes from
/// `Region::signed_distance_grad`.
///
/// # Example
/// ```
/// use molpack::region::{InsideSphereRegion, RegionExt};
/// use molpack::{RegionRestraint, AtomRestraint};
///
/// let shell = InsideSphereRegion::new([0.0; 3], 10.0)
///     .and(InsideSphereRegion::new([0.0; 3], 5.0).not());
/// let restraint = RegionRestraint(shell);
/// // At x=(7,0,0) the shell is satisfied → f == 0
/// assert_eq!(restraint.f(&[7.0, 0.0, 0.0], 1.0, 1.0), 0.0);
/// // At x=(15,0,0) we are outside the outer sphere → f > 0
/// assert!(restraint.f(&[15.0, 0.0, 0.0], 1.0, 1.0) > 0.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct RegionRestraint<R: Region>(pub R);

impl<R: Region + 'static> AtomRestraint for RegionRestraint<R> {
    fn f(&self, x: &[F; 3], _scale: F, scale2: F) -> F {
        let d = self.0.signed_distance(x);
        let v = d.max(0.0);
        scale2 * v * v
    }

    fn fg(&self, x: &[F; 3], scale: F, scale2: F, g: &mut [F; 3]) -> F {
        let d = self.0.signed_distance(x);
        if d > 0.0 {
            let grad = self.0.signed_distance_grad(x);
            let coeff = 2.0 * scale2 * d;
            g[0] += coeff * grad[0];
            g[1] += coeff * grad[1];
            g[2] += coeff * grad[2];
        }
        self.f(x, scale, scale2)
    }

    fn declared_cell(&self) -> Option<CellDeclaration> {
        self.0.declared_cell()
    }

    fn is_closed_mesh(&self) -> bool {
        self.0.is_closed_mesh()
    }
}

// ============================================================================
// Concrete regions (starting set — more can be added incrementally)
// ============================================================================

/// Axis-aligned box region (inside test).
#[derive(Debug, Clone, Copy)]
pub struct InsideBoxRegion {
    pub min: [F; 3],
    pub max: [F; 3],
}

impl InsideBoxRegion {
    pub fn new(min: [F; 3], max: [F; 3]) -> Self {
        Self { min, max }
    }
}

impl Region for InsideBoxRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        (0..3).all(|k| x[k] >= self.min[k] && x[k] <= self.max[k])
    }

    fn signed_distance(&self, x: &[F; 3]) -> F {
        // Signed distance to an axis-aligned box:
        // d = max_k max(min_k - x_k, x_k - max_k)
        // Negative inside (both terms negative), positive outside.
        let mut d = F::NEG_INFINITY;
        for ((xk, &lo_k), &hi_k) in x.iter().zip(self.min.iter()).zip(self.max.iter()) {
            d = d.max(lo_k - xk).max(xk - hi_k);
        }
        d
    }

    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        // The max is attained on one face; gradient points outward from that face.
        let mut best_d = F::NEG_INFINITY;
        let mut best_axis = 0usize;
        let mut best_sign = 0.0 as F;
        for (k, ((xk, &lo_k), &hi_k)) in x
            .iter()
            .zip(self.min.iter())
            .zip(self.max.iter())
            .enumerate()
        {
            let lo = lo_k - xk; // gradient component -1 on axis k
            let hi = xk - hi_k; // gradient component +1 on axis k
            if lo > best_d {
                best_d = lo;
                best_axis = k;
                best_sign = -1.0;
            }
            if hi > best_d {
                best_d = hi;
                best_axis = k;
                best_sign = 1.0;
            }
        }
        let mut g = [0.0; 3];
        g[best_axis] = best_sign;
        g
    }

    fn bounding_box(&self) -> Option<Aabb> {
        Some(Aabb {
            min: self.min,
            max: self.max,
        })
    }
}

/// Spherical region (inside test).
#[derive(Debug, Clone, Copy)]
pub struct InsideSphereRegion {
    pub center: [F; 3],
    pub radius: F,
}

impl InsideSphereRegion {
    pub fn new(center: [F; 3], radius: F) -> Self {
        Self { center, radius }
    }
}

impl Region for InsideSphereRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        let c = self.center;
        let d2 = (x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2);
        d2 <= self.radius.powi(2)
    }

    fn signed_distance(&self, x: &[F; 3]) -> F {
        let c = self.center;
        let d = ((x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2)).sqrt();
        d - self.radius
    }

    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let c = self.center;
        let (dx, dy, dz) = (x[0] - c[0], x[1] - c[1], x[2] - c[2]);
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < 1e-12 {
            [0.0; 3]
        } else {
            [dx / d, dy / d, dz / d]
        }
    }

    fn bounding_box(&self) -> Option<Aabb> {
        let r = self.radius;
        let c = self.center;
        Some(Aabb {
            min: [c[0] - r, c[1] - r, c[2] - r],
            max: [c[0] + r, c[1] + r, c[2] + r],
        })
    }
}

/// Outside-sphere region: `!InsideSphereRegion` with a bespoke impl so
/// `signed_distance` avoids the `Not(...)` double negation at call time.
#[derive(Debug, Clone, Copy)]
pub struct OutsideSphereRegion {
    pub center: [F; 3],
    pub radius: F,
}

impl OutsideSphereRegion {
    pub fn new(center: [F; 3], radius: F) -> Self {
        Self { center, radius }
    }
}

impl Region for OutsideSphereRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        let c = self.center;
        let d2 = (x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2);
        d2 >= self.radius.powi(2)
    }

    fn signed_distance(&self, x: &[F; 3]) -> F {
        let c = self.center;
        let d = ((x[0] - c[0]).powi(2) + (x[1] - c[1]).powi(2) + (x[2] - c[2]).powi(2)).sqrt();
        self.radius - d
    }

    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let c = self.center;
        let (dx, dy, dz) = (x[0] - c[0], x[1] - c[1], x[2] - c[2]);
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < 1e-12 {
            [0.0; 3]
        } else {
            [-dx / d, -dy / d, -dz / d]
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

// ============================================================================
// Primitive cell
// ============================================================================

/// The primitive cell of a lattice, or a fractional sub-slice of it.
///
/// This is the region an orthorhombic packer cannot express. Membership is
/// decided in fractional coordinates, so the cell may be hexagonal, monoclinic
/// or fully triclinic; distances are reported in Ångström, measured
/// perpendicular to the bounding lattice planes, so the penalty is comparable
/// with every other restraint regardless of how tilted the cell is.
///
/// Declaring this region also declares the lattice — see
/// [`Region::declared_cell`] — so a caller writes the cell once.
#[derive(Debug, Clone)]
pub struct InsideCellRegion {
    bx: SimBox,
    lo: [F; 3],
    hi: [F; 3],
    /// Interplanar spacing along each reciprocal direction: the factor turning
    /// a fractional offset into an Ångström distance from that pair of faces.
    spacing: [F; 3],
    /// Unit outward normal of each face pair, i.e. the normalised rows of H⁻¹.
    normal: [[F; 3]; 3],
    pbc: [bool; 3],
}

impl InsideCellRegion {
    /// Cell from lengths (Å) and angles (degrees).
    pub fn from_lengths_angles(
        lengths: [F; 3],
        angles_deg: [F; 3],
        pbc: [bool; 3],
    ) -> Option<Self> {
        let h = SimBox::matrix_from_lengths_angles(lengths, angles_deg).ok()?;
        Self::from_simbox(SimBox::new(h, array![0.0, 0.0, 0.0], pbc).ok()?)
    }

    /// Cell from a lattice matrix whose **columns** are the lattice vectors.
    pub fn from_matrix(h: [[F; 3]; 3], origin: [F; 3], pbc: [bool; 3]) -> Option<Self> {
        let matrix = array![
            [h[0][0], h[0][1], h[0][2]],
            [h[1][0], h[1][1], h[1][2]],
            [h[2][0], h[2][1], h[2][2]]
        ];
        Self::from_simbox(SimBox::new(matrix, array![origin[0], origin[1], origin[2]], pbc).ok()?)
    }

    /// Cell from an existing [`SimBox`].
    pub fn from_simbox(bx: SimBox) -> Option<Self> {
        let inv = bx.inv_view();
        let mut spacing = [0.0; 3];
        let mut normal = [[0.0; 3]; 3];
        for k in 0..3 {
            let row = [inv[[k, 0]], inv[[k, 1]], inv[[k, 2]]];
            let norm = (row[0] * row[0] + row[1] * row[1] + row[2] * row[2]).sqrt();
            if norm <= 0.0 || !norm.is_finite() {
                return None;
            }
            spacing[k] = 1.0 / norm;
            normal[k] = [row[0] / norm, row[1] / norm, row[2] / norm];
        }
        let pbc = bx.pbc();
        Some(Self {
            bx,
            lo: [0.0; 3],
            hi: [1.0; 3],
            spacing,
            normal,
            pbc,
        })
    }

    /// Restrict to a fractional sub-slice, e.g. the middle third along `c`.
    ///
    /// Bounds outside `[0, 1]` are accepted — a slab may deliberately sit
    /// proud of the cell — but `lo` must stay below `hi` on every axis.
    pub fn with_fractional_bounds(mut self, lo: [F; 3], hi: [F; 3]) -> Option<Self> {
        if (0..3).any(|k| lo[k] >= hi[k] || !lo[k].is_finite() || !hi[k].is_finite()) {
            return None;
        }
        self.lo = lo;
        self.hi = hi;
        Some(self)
    }

    /// Signed distance to the nearest bounding plane, and which face won.
    ///
    /// Returns `(distance, axis, sign)` where `sign` is `+1` when the point is
    /// past the `hi` face and `-1` when past the `lo` face.
    ///
    /// Every axis carries faces, including periodic ones. Under periodicity a
    /// molecule at fractional `1.04` is the same configuration as one at
    /// `0.04` — the pair kernel's minimum image cannot tell them apart — so the
    /// excursion is not a physical defect. It is still confined, because a
    /// packer has to emit coordinates someone can use: with nothing holding
    /// them, molecules drift across hundreds of lattice images (measured: ±1500
    /// Å from a ±`sidemax` initial placement) and every consumer then has to
    /// wrap before the result means anything.
    ///
    /// Because the penalty is quadratic rather than a hard wall, equilibrium
    /// leaves sub-tolerance excursions past a face — the same behaviour as
    /// every other `Inside*` restraint.
    #[inline]
    fn nearest_face(&self, x: &[F; 3]) -> (F, usize, F) {
        let f = self.bx.make_fractional_raw_arr3(*x);
        let mut best = (F::NEG_INFINITY, 0usize, 1.0 as F);
        for (k, &fk) in f.iter().enumerate() {
            let below = (self.lo[k] - fk) * self.spacing[k];
            if below > best.0 {
                best = (below, k, -1.0);
            }
            let above = (fk - self.hi[k]) * self.spacing[k];
            if above > best.0 {
                best = (above, k, 1.0);
            }
        }
        best
    }
}

impl Region for InsideCellRegion {
    fn contains(&self, x: &[F; 3]) -> bool {
        self.nearest_face(x).0 <= 0.0
    }

    fn signed_distance(&self, x: &[F; 3]) -> F {
        self.nearest_face(x).0
    }

    fn signed_distance_grad(&self, x: &[F; 3]) -> [F; 3] {
        let (_, axis, sign) = self.nearest_face(x);
        let n = self.normal[axis];
        [sign * n[0], sign * n[1], sign * n[2]]
    }

    fn declared_cell(&self) -> Option<CellDeclaration> {
        let h = self.bx.h_view();
        let o = self.bx.origin_view();
        Some((
            [
                [h[[0, 0]], h[[0, 1]], h[[0, 2]]],
                [h[[1, 0]], h[[1, 1]], h[[1, 2]]],
                [h[[2, 0]], h[[2, 1]], h[[2, 2]]],
            ],
            [o[0], o[1], o[2]],
            self.pbc,
        ))
    }

    fn bounding_box(&self) -> Option<Aabb> {
        let mut min = [F::INFINITY; 3];
        let mut max = [F::NEG_INFINITY; 3];
        for i in 0..8 {
            let frac = array![
                if i & 1 == 0 { self.lo[0] } else { self.hi[0] },
                if i & 2 == 0 { self.lo[1] } else { self.hi[1] },
                if i & 4 == 0 { self.lo[2] } else { self.hi[2] },
            ];
            let corner = self.bx.make_cartesian(frac.view());
            for k in 0..3 {
                min[k] = min[k].min(corner[k]);
                max[k] = max[k].max(corner[k]);
            }
        }
        Some(Aabb { min, max })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: F = 1e-6;

    // ── boolean algebra laws ────────────────────────────────────────────────

    #[test]
    fn and_is_intersection() {
        let a = InsideBoxRegion::new([0.0; 3], [10.0; 3]);
        let b = InsideSphereRegion::new([5.0; 3], 4.0);
        let c = And(a, b);
        // Inside both
        assert!(c.contains(&[5.0, 5.0, 5.0]));
        // Inside box only
        assert!(!c.contains(&[1.0, 1.0, 1.0]));
        // Inside sphere only — impossible here since sphere ⊂ box
        // Outside both
        assert!(!c.contains(&[-1.0, -1.0, -1.0]));
    }

    #[test]
    fn or_is_union() {
        let a = InsideSphereRegion::new([0.0; 3], 5.0);
        let b = InsideSphereRegion::new([20.0, 0.0, 0.0], 5.0);
        let u = Or(a, b);
        assert!(u.contains(&[0.0, 0.0, 0.0]));
        assert!(u.contains(&[20.0, 0.0, 0.0]));
        assert!(!u.contains(&[10.0, 0.0, 0.0]));
    }

    #[test]
    fn not_is_complement() {
        let a = InsideSphereRegion::new([0.0; 3], 5.0);
        let n = Not(a);
        assert!(!n.contains(&[0.0, 0.0, 0.0]));
        assert!(n.contains(&[10.0, 0.0, 0.0]));
    }

    #[test]
    fn de_morgan_and() {
        // !(A ∧ B) == !A ∨ !B
        let a = InsideSphereRegion::new([0.0; 3], 5.0);
        let b = InsideBoxRegion::new([-3.0; 3], [3.0; 3]);
        let lhs = Not(And(a, b));
        let rhs = Or(Not(a), Not(b));
        for pt in &[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        ] {
            assert_eq!(
                lhs.contains(pt),
                rhs.contains(pt),
                "de Morgan mismatch at {pt:?}"
            );
        }
    }

    // ── signed_distance sign correctness ─────────────────────────────────────

    #[test]
    fn sphere_signed_distance_sign() {
        let s = InsideSphereRegion::new([0.0; 3], 5.0);
        assert!(s.signed_distance(&[0.0, 0.0, 0.0]) < 0.0);
        assert!((s.signed_distance(&[5.0, 0.0, 0.0])).abs() < TOL);
        assert!(s.signed_distance(&[10.0, 0.0, 0.0]) > 0.0);
    }

    #[test]
    fn box_signed_distance_sign() {
        let b = InsideBoxRegion::new([0.0; 3], [10.0; 3]);
        assert!(b.signed_distance(&[5.0, 5.0, 5.0]) < 0.0);
        assert!(b.signed_distance(&[15.0, 5.0, 5.0]) > 0.0);
        assert!(b.signed_distance(&[-5.0, 5.0, 5.0]) > 0.0);
    }

    #[test]
    fn contains_matches_signed_distance() {
        let regions: Vec<Box<dyn Region>> = vec![
            Box::new(InsideSphereRegion::new([0.0; 3], 3.0)),
            Box::new(InsideBoxRegion::new([0.0; 3], [5.0; 3])),
            Box::new(OutsideSphereRegion::new([0.0; 3], 2.0)),
        ];
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [3.0, 3.0, 3.0],
            [-1.0, 0.0, 0.0],
        ];
        for r in &regions {
            for pt in &pts {
                let sd = r.signed_distance(pt);
                // contains ⇔ sd <= 0 (within tolerance for boundary)
                if sd < -TOL {
                    assert!(r.contains(pt), "sd={sd}<0 but !contains at {pt:?}");
                }
                if sd > TOL {
                    assert!(!r.contains(pt), "sd={sd}>0 but contains at {pt:?}");
                }
            }
        }
    }

    // ── RegionRestraint lifts to AtomRestraint ───────────────────────────────────

    #[test]
    fn from_region_penalty_zero_inside() {
        let r = RegionRestraint(InsideSphereRegion::new([0.0; 3], 5.0));
        assert_eq!(r.f(&[0.0, 0.0, 0.0], 1.0, 1.0), 0.0);
    }

    #[test]
    fn from_region_penalty_positive_outside() {
        let r = RegionRestraint(InsideSphereRegion::new([0.0; 3], 5.0));
        assert!(r.f(&[10.0, 0.0, 0.0], 1.0, 1.0) > 0.0);
    }

    #[test]
    fn from_region_gradient_matches_finite_diff() {
        let r = RegionRestraint(
            InsideBoxRegion::new([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
                .and(Not(InsideSphereRegion::new([5.0, 5.0, 5.0], 2.0))),
        );
        // Point outside the box in x → composite region reports "outside" via box face
        let x = [15.0, 5.0, 5.0];
        let mut g = [0.0; 3];
        let _ = r.fg(&x, 1.0, 1.0, &mut g);

        // central finite difference
        let h: F = 1e-5;
        for k in 0..3 {
            let mut xp = x;
            xp[k] += h;
            let mut xm = x;
            xm[k] -= h;
            let fd = (r.f(&xp, 1.0, 1.0) - r.f(&xm, 1.0, 1.0)) / (2.0 * h);
            assert!(
                (g[k] - fd).abs() < 1e-3,
                "gradient mismatch axis {k}: analytic={} fd={} err={}",
                g[k],
                fd,
                (g[k] - fd).abs()
            );
        }
    }

    #[test]
    fn region_ext_chain() {
        // RegionExt provides ergonomic .and() / .not()
        let shell = InsideSphereRegion::new([0.0; 3], 10.0)
            .and(InsideSphereRegion::new([0.0; 3], 5.0).not());
        assert!(shell.contains(&[7.0, 0.0, 0.0]));
        assert!(!shell.contains(&[3.0, 0.0, 0.0]));
        assert!(!shell.contains(&[15.0, 0.0, 0.0]));
    }

    #[test]
    fn gradient_accumulates_not_overwrite() {
        let r = RegionRestraint(InsideSphereRegion::new([0.0; 3], 1.0));
        let mut g = [100.0; 3];
        let _ = r.fg(&[2.0, 0.0, 0.0], 1.0, 1.0, &mut g);
        assert!(g[0] > 100.0, "should accumulate");
        assert!((g[1] - 100.0).abs() < TOL);
        assert!((g[2] - 100.0).abs() < TOL);
    }
}
