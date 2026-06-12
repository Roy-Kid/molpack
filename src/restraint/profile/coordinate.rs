//! Reaction coordinate ξ(x) for the profile-restraint family.
//!
//! A [`Coordinate`] maps a Cartesian position `x` (Å) to a scalar reaction
//! coordinate ξ (a length in Å) and supplies the analytic gradient ∇ξ. A later
//! composing restraint applies the chain rule `∇ₓU = (dU/dξ)·∇ξ`; this module is
//! pure geometry — no penalty/distribution math and no `Restraint` impl.
//!
//! Four flavours are provided, each minimum-image aware (the reference
//! plane/centre/axis is FIXED in the box, so the `x − reference` delta is wrapped
//! per active periodic axis before any projection/norm) and each guarded so the
//! gradient stays finite even at its geometric singularity:
//!
//! - **planar** — ξ = n̂·(x−x₀), ∇ξ = n̂ (constant; n̂ unit on construction; no
//!   singularity).
//! - **radial** — ξ = ‖x−c‖, ∇ξ = (x−c)/‖x−c‖ = r̂ (undefined at the centre;
//!   inner clamp `r_guard`).
//! - **cylindrical** — ρ⃗ = (x−a) − ((x−a)·û)û, ξ = ‖ρ⃗‖, ∇ξ = ρ⃗/‖ρ⃗‖. The
//!   projection `(I−ûûᵀ)` is symmetric idempotent, so ∇ξ = ρ⃗/‖ρ⃗‖ (û unit on
//!   construction; undefined on the axis; inner clamp `r_guard`).
//! - **region-distance** — ξ = `Region::signed_distance(x)`, ∇ξ =
//!   `Region::signed_distance_grad(x)` (outward unit normal a.e.; ignores PBC,
//!   relying on the region's own convention).

use std::sync::Arc;

use molrs::types::F;

use crate::cell::delta_vector;
use crate::region::Region;

/// Minimum-image carrier passed by the caller. Mirrors the PBC representation of
/// [`delta_vector`](crate::cell::delta_vector): box edge lengths plus a per-axis
/// periodic flag. A non-periodic axis (`periodic[k] == false`) is never wrapped,
/// so `PbcWrap::none()` reproduces raw Cartesian differences.
#[derive(Debug, Clone, Copy)]
pub struct PbcWrap {
    /// Box edge length along each axis (Å). Ignored on non-periodic axes.
    pub pbc_length: [F; 3],
    /// Per-axis periodic flag; only `true` axes are minimum-image wrapped.
    pub periodic: [bool; 3],
}

impl PbcWrap {
    /// A fully non-periodic carrier: every `x − reference` delta is the raw
    /// Cartesian difference.
    pub fn none() -> Self {
        Self {
            pbc_length: [0.0; 3],
            periodic: [false; 3],
        }
    }

    /// A carrier with the given box lengths and per-axis periodic flags.
    pub fn new(pbc_length: [F; 3], periodic: [bool; 3]) -> Self {
        Self {
            pbc_length,
            periodic,
        }
    }

    /// Minimum-image `x − reference` honoring the active periodic axes.
    #[inline]
    fn delta(&self, x: &[F; 3], reference: &[F; 3]) -> [F; 3] {
        delta_vector(x, reference, &self.pbc_length, &self.periodic)
    }
}

/// Error returned by fallible [`Coordinate`] constructors when a direction
/// parameter is degenerate (near-zero length) and cannot be normalized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateError {
    /// The planar normal had near-zero length and could not be normalized.
    ZeroNormal,
    /// The cylindrical axis direction had near-zero length and could not be
    /// normalized.
    ZeroAxis,
}

impl std::fmt::Display for CoordinateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoordinateError::ZeroNormal => {
                write!(f, "planar coordinate normal has near-zero length")
            }
            CoordinateError::ZeroAxis => {
                write!(f, "cylindrical coordinate axis has near-zero length")
            }
        }
    }
}

impl std::error::Error for CoordinateError {}

/// Below this length a direction vector is treated as degenerate and rejected at
/// construction (planar normal / cylindrical axis).
const MIN_DIR_LEN: F = 1e-12;

/// Reaction coordinate ξ(x) for the profile-restraint family.
///
/// Immutable: constructed once from validated parameters, holds no mutable
/// state, and is shared read-only across packing iterations (`Send + Sync`).
/// `Clone` is available because every payload is `Clone` (the region variant
/// holds an `Arc`).
#[derive(Debug, Clone)]
pub enum Coordinate {
    /// Signed distance to a plane: ξ = n̂·(x−x₀), ∇ξ = n̂. `normal` is unit on
    /// construction, so there is no singularity.
    Planar {
        /// Unit plane normal n̂.
        normal: [F; 3],
        /// A point x₀ on the plane (Å).
        point: [F; 3],
    },
    /// Distance from a centre: ξ = ‖x−c‖, ∇ξ = (x−c)/‖x−c‖. Undefined at the
    /// centre; inside `r_guard` ξ is reported as `r_guard` and ∇ξ frozen to a
    /// finite unit vector.
    Radial {
        /// Centre c (Å).
        center: [F; 3],
        /// Inner clamp radius (Å): for r < `r_guard`, ξ = `r_guard`.
        r_guard: F,
    },
    /// Perpendicular distance to an axis: ρ⃗ = (x−a) − ((x−a)·û)û, ξ = ‖ρ⃗‖,
    /// ∇ξ = ρ⃗/‖ρ⃗‖. `axis_dir` is unit on construction; undefined on the axis;
    /// inside `r_guard` ξ is reported as `r_guard` and ∇ξ frozen.
    Cylindrical {
        /// A point a on the axis (Å).
        axis_origin: [F; 3],
        /// Unit axis direction û.
        axis_dir: [F; 3],
        /// Inner clamp radius (Å): for s < `r_guard`, ξ = `r_guard`.
        r_guard: F,
    },
    /// Signed distance to a region surface: ξ = `signed_distance(x)`,
    /// ∇ξ = `signed_distance_grad(x)`. Ignores PBC, relying on the region's own
    /// convention.
    RegionDistance(Arc<dyn Region>),
}

impl Coordinate {
    /// Construct a planar coordinate, normalizing `normal` to unit length.
    /// Returns [`CoordinateError::ZeroNormal`] when `normal` is degenerate.
    pub fn planar(normal: [F; 3], point: [F; 3]) -> Result<Self, CoordinateError> {
        let normal = normalize(&normal).ok_or(CoordinateError::ZeroNormal)?;
        Ok(Coordinate::Planar { normal, point })
    }

    /// Construct a radial coordinate. Infallible. `r_guard` is the inner clamp
    /// radius (its absolute value is used; negative inputs are folded).
    pub fn radial(center: [F; 3], r_guard: F) -> Self {
        Coordinate::Radial {
            center,
            r_guard: r_guard.abs(),
        }
    }

    /// Construct a cylindrical coordinate, normalizing `axis_dir` to unit length.
    /// Returns [`CoordinateError::ZeroAxis`] when `axis_dir` is degenerate.
    pub fn cylindrical(
        axis_origin: [F; 3],
        axis_dir: [F; 3],
        r_guard: F,
    ) -> Result<Self, CoordinateError> {
        let axis_dir = normalize(&axis_dir).ok_or(CoordinateError::ZeroAxis)?;
        Ok(Coordinate::Cylindrical {
            axis_origin,
            axis_dir,
            r_guard: r_guard.abs(),
        })
    }

    /// Construct a region-distance coordinate. Infallible.
    pub fn region_distance(region: Arc<dyn Region>) -> Self {
        Coordinate::RegionDistance(region)
    }

    /// Evaluate the reaction coordinate ξ (Å) at `x`.
    ///
    /// For planar/radial/cylindrical the `x − reference` delta is minimum-image
    /// wrapped via `pbc` before projecting/norming; region-distance ignores
    /// `pbc`.
    pub fn xi(&self, x: &[F; 3], pbc: &PbcWrap) -> F {
        match self {
            Coordinate::Planar { normal, point } => {
                let d = pbc.delta(x, point);
                dot(normal, &d)
            }
            Coordinate::Radial { center, r_guard } => {
                let d = pbc.delta(x, center);
                norm(&d).max(*r_guard)
            }
            Coordinate::Cylindrical {
                axis_origin,
                axis_dir,
                r_guard,
            } => {
                let rho = perp_component(&pbc.delta(x, axis_origin), axis_dir);
                norm(&rho).max(*r_guard)
            }
            Coordinate::RegionDistance(region) => region.signed_distance(x),
        }
    }

    /// Evaluate the analytic gradient ∇ξ at `x`.
    ///
    /// Always finite: the radial/cylindrical inner clamp freezes ∇ξ to a finite
    /// unit vector (a fixed fallback direction exactly on the singularity), so
    /// the returned gradient norm is bounded and never NaN/Inf.
    pub fn grad_xi(&self, x: &[F; 3], pbc: &PbcWrap) -> [F; 3] {
        match self {
            Coordinate::Planar { normal, .. } => *normal,
            Coordinate::Radial { center, r_guard } => {
                let d = pbc.delta(x, center);
                unit_or_fallback(&d, norm(&d), *r_guard)
            }
            Coordinate::Cylindrical {
                axis_origin,
                axis_dir,
                r_guard,
            } => {
                let rho = perp_component(&pbc.delta(x, axis_origin), axis_dir);
                // Fallback must lie in the plane ⟂ axis; FALLBACK_DIR may be
                // parallel to the axis, so project it before use.
                unit_or_fallback_perp(&rho, axis_dir, *r_guard)
            }
            Coordinate::RegionDistance(region) => region.signed_distance_grad(x),
        }
    }
}

/// Fixed fallback unit direction used exactly on a radial/cylindrical
/// singularity, where r̂ is undefined. Any unit vector keeps ∇ξ finite; the
/// magnitude (1) is what the clamp guarantees.
const FALLBACK_DIR: [F; 3] = [1.0, 0.0, 0.0];

#[inline]
fn dot(a: &[F; 3], b: &[F; 3]) -> F {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(v: &[F; 3]) -> F {
    dot(v, v).sqrt()
}

/// Normalize `v`; returns `None` if `v` is shorter than [`MIN_DIR_LEN`].
#[inline]
fn normalize(v: &[F; 3]) -> Option<[F; 3]> {
    let n = norm(v);
    if n < MIN_DIR_LEN {
        None
    } else {
        Some([v[0] / n, v[1] / n, v[2] / n])
    }
}

/// Component of `d` perpendicular to the unit axis `u`: ρ⃗ = d − (d·û)û.
#[inline]
fn perp_component(d: &[F; 3], u: &[F; 3]) -> [F; 3] {
    let proj = dot(d, u);
    [d[0] - proj * u[0], d[1] - proj * u[1], d[2] - proj * u[2]]
}

/// Unit vector `v/‖v‖` when `len >= r_guard`; otherwise the fixed fallback unit
/// direction. Keeps ∇ξ finite and unit-norm at the singularity.
#[inline]
fn unit_or_fallback(v: &[F; 3], len: F, r_guard: F) -> [F; 3] {
    if len >= r_guard && len >= MIN_DIR_LEN {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        FALLBACK_DIR
    }
}

/// Like [`unit_or_fallback`] but for the cylindrical case: the fallback must lie
/// in the plane perpendicular to the axis `u`. If the fixed fallback is parallel
/// to the axis, a perpendicular substitute is chosen.
#[inline]
fn unit_or_fallback_perp(rho: &[F; 3], u: &[F; 3], r_guard: F) -> [F; 3] {
    let len = norm(rho);
    if len >= r_guard && len >= MIN_DIR_LEN {
        return [rho[0] / len, rho[1] / len, rho[2] / len];
    }
    // On the axis: pick any unit vector ⟂ û. Project the fixed fallback onto the
    // perpendicular plane; if it is (near-)parallel to û, use a different basis
    // vector that cannot also be parallel.
    let candidate = perp_component(&FALLBACK_DIR, u);
    if let Some(unit) = normalize(&candidate) {
        unit
    } else {
        // FALLBACK_DIR was parallel to the axis; [0,1,0] cannot also be.
        normalize(&perp_component(&[0.0, 1.0, 0.0], u)).unwrap_or(FALLBACK_DIR)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::region::InsideBoxRegion;

    const TOL: F = 1e-6;

    /// 3-point central finite difference of `xi` along each axis.
    fn fd_grad(c: &Coordinate, x: &[F; 3], pbc: &PbcWrap, h: F) -> [F; 3] {
        let mut g = [0.0; 3];
        for k in 0..3 {
            let mut xp = *x;
            xp[k] += h;
            let mut xm = *x;
            xm[k] -= h;
            g[k] = (c.xi(&xp, pbc) - c.xi(&xm, pbc)) / (2.0 * h);
        }
        g
    }

    /// ac-001: analytic ∇ξ matches central FD componentwise for every flavour.
    fn assert_grad_matches_fd(c: &Coordinate, x: &[F; 3], pbc: &PbcWrap, label: &str) {
        let analytic = c.grad_xi(x, pbc);
        let fd = fd_grad(c, x, pbc, 1e-6);
        for k in 0..3 {
            assert!(
                (analytic[k] - fd[k]).abs() <= TOL,
                "{label}: ∇ξ axis {k} analytic={} fd={} err={}",
                analytic[k],
                fd[k],
                (analytic[k] - fd[k]).abs()
            );
        }
    }

    // ── ac-001: analytic vs finite-difference for all four flavours ──────────

    #[test]
    fn planar_grad_matches_fd() {
        let c = Coordinate::planar([1.0, 2.0, -3.0], [0.5, -1.0, 2.0]).unwrap();
        assert_grad_matches_fd(&c, &[3.0, 4.0, 5.0], &PbcWrap::none(), "planar");
    }

    #[test]
    fn planar_xi_is_closed_form() {
        // Unit normal along +x, plane through origin → ξ = x-coordinate.
        let c = Coordinate::planar([2.0, 0.0, 0.0], [0.0, 0.0, 0.0]).unwrap();
        let x = [3.0, 9.0, -4.0];
        assert!((c.xi(&x, &PbcWrap::none()) - 3.0).abs() < TOL);
        // ∇ξ is the constant unit normal.
        assert_eq!(c.grad_xi(&x, &PbcWrap::none()), [1.0, 0.0, 0.0]);
    }

    #[test]
    fn radial_grad_matches_fd() {
        let c = Coordinate::radial([1.0, -2.0, 0.5], 0.1);
        assert_grad_matches_fd(&c, &[4.0, 1.0, 3.0], &PbcWrap::none(), "radial");
    }

    #[test]
    fn radial_xi_is_closed_form() {
        let c = Coordinate::radial([0.0, 0.0, 0.0], 0.1);
        let x = [3.0, 4.0, 0.0];
        assert!((c.xi(&x, &PbcWrap::none()) - 5.0).abs() < TOL);
        // ∇ξ unit norm off-singularity.
        let g = c.grad_xi(&x, &PbcWrap::none());
        assert!((norm(&g) - 1.0).abs() < TOL);
    }

    #[test]
    fn cylindrical_grad_matches_fd() {
        let c = Coordinate::cylindrical([0.0, 0.0, 0.0], [0.0, 0.0, 2.0], 0.1).unwrap();
        // Off the axis (axis is z): perpendicular distance is in the xy-plane.
        assert_grad_matches_fd(&c, &[3.0, 4.0, 7.0], &PbcWrap::none(), "cylindrical");
    }

    #[test]
    fn cylindrical_xi_is_closed_form() {
        // Axis = z through origin → ξ = sqrt(x²+y²), independent of z.
        let c = Coordinate::cylindrical([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.1).unwrap();
        let x = [3.0, 4.0, 100.0];
        assert!((c.xi(&x, &PbcWrap::none()) - 5.0).abs() < TOL);
        let g = c.grad_xi(&x, &PbcWrap::none());
        assert!((norm(&g) - 1.0).abs() < TOL);
        // Gradient lies in the xy-plane (⟂ axis).
        assert!(g[2].abs() < TOL);
    }

    #[test]
    fn region_distance_grad_matches_fd() {
        let region = Arc::new(InsideBoxRegion::new([0.0; 3], [10.0; 3]));
        let c = Coordinate::region_distance(region);
        // Outside the box on +x: nearest face is x=10, signed distance smooth here.
        assert_grad_matches_fd(&c, &[13.0, 5.0, 5.0], &PbcWrap::none(), "region");
    }

    #[test]
    fn region_distance_xi_is_signed_distance() {
        let region = Arc::new(InsideBoxRegion::new([0.0; 3], [10.0; 3]));
        let c = Coordinate::region_distance(region.clone());
        let x = [13.0, 5.0, 5.0];
        assert!((c.xi(&x, &PbcWrap::none()) - region.signed_distance(&x)).abs() < TOL);
        // ‖∇d‖ = 1 a.e. (away from medial axis).
        let g = c.grad_xi(&x, &PbcWrap::none());
        assert!((norm(&g) - 1.0).abs() < TOL);
    }

    // ── ac-002: finite, bounded gradient inside r_guard ──────────────────────

    #[test]
    fn radial_at_centre_is_finite() {
        let r_guard = 0.25;
        let c = Coordinate::radial([1.0, 2.0, 3.0], r_guard);
        let centre = [1.0, 2.0, 3.0];
        // ξ clamped to r_guard.
        assert!((c.xi(&centre, &PbcWrap::none()) - r_guard).abs() < TOL);
        let g = c.grad_xi(&centre, &PbcWrap::none());
        for (k, gk) in g.iter().enumerate() {
            assert!(gk.is_finite(), "radial ∇ξ axis {k} not finite: {gk}");
        }
        assert!(
            norm(&g) <= 1.0 + TOL,
            "radial ∇ξ norm unbounded: {}",
            norm(&g)
        );
    }

    #[test]
    fn cylindrical_on_axis_is_finite() {
        let r_guard = 0.25;
        let c = Coordinate::cylindrical([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], r_guard).unwrap();
        // A point exactly on the z-axis.
        let on_axis = [0.0, 0.0, 7.0];
        assert!((c.xi(&on_axis, &PbcWrap::none()) - r_guard).abs() < TOL);
        let g = c.grad_xi(&on_axis, &PbcWrap::none());
        for (k, gk) in g.iter().enumerate() {
            assert!(gk.is_finite(), "cyl ∇ξ axis {k} not finite: {gk}");
        }
        assert!(norm(&g) <= 1.0 + TOL, "cyl ∇ξ norm unbounded: {}", norm(&g));
        // Fallback must be perpendicular to the axis (z).
        assert!(g[2].abs() < TOL, "cyl fallback not ⟂ axis: {:?}", g);
    }

    #[test]
    fn radial_just_inside_r_guard_reports_r_guard() {
        let r_guard = 0.5;
        let c = Coordinate::radial([0.0, 0.0, 0.0], r_guard);
        // Radius 0.3 < r_guard → ξ reported as r_guard.
        let x = [0.3, 0.0, 0.0];
        assert!((c.xi(&x, &PbcWrap::none()) - r_guard).abs() < TOL);
    }

    // ── ac-003: minimum-image wrap for planar/radial/cylindrical ─────────────

    #[test]
    fn radial_uses_minimum_image() {
        // Box length 10 on x (periodic), reference at x=9, query at x=1.
        // Raw delta = 1-9 = -8 → ‖8‖. Wrapped: -8 + 10 = 2 → ‖2‖.
        let pbc = PbcWrap::new([10.0, 10.0, 10.0], [true, false, false]);
        let c = Coordinate::radial([9.0, 0.0, 0.0], 0.0);
        let x = [1.0, 0.0, 0.0];
        let xi = c.xi(&x, &pbc);
        assert!(
            (xi - 2.0).abs() < TOL,
            "wrapped radial ξ={xi}, expected 2.0"
        );
        // Without PBC it would be the raw 8.0.
        let xi_raw = c.xi(&x, &PbcWrap::none());
        assert!(
            (xi_raw - 8.0).abs() < TOL,
            "raw radial ξ={xi_raw}, expected 8.0"
        );
    }

    #[test]
    fn planar_uses_minimum_image() {
        // Normal +x, plane at x=9, periodic box length 10. Query at x=1:
        // raw n̂·(x−x₀) = 1−9 = −8; wrapped = +2.
        let pbc = PbcWrap::new([10.0, 10.0, 10.0], [true, false, false]);
        let c = Coordinate::planar([1.0, 0.0, 0.0], [9.0, 0.0, 0.0]).unwrap();
        let x = [1.0, 0.0, 0.0];
        assert!((c.xi(&x, &pbc) - 2.0).abs() < TOL, "wrapped planar ξ wrong");
        assert!(
            (c.xi(&x, &PbcWrap::none()) - (-8.0)).abs() < TOL,
            "raw planar ξ wrong"
        );
    }

    #[test]
    fn cylindrical_uses_minimum_image() {
        // Axis = z through (9,0,*), periodic on x with length 10. Query (1,3,5):
        // raw delta x = 1−9 = −8 → ρ in xy = (−8,3) → ‖·‖ = sqrt(73).
        // wrapped delta x = +2 → ρ = (2,3) → ‖·‖ = sqrt(13).
        let pbc = PbcWrap::new([10.0, 10.0, 10.0], [true, false, false]);
        let c = Coordinate::cylindrical([9.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0).unwrap();
        let x = [1.0, 3.0, 5.0];
        let expected = (2.0_f64 * 2.0 + 3.0 * 3.0).sqrt();
        assert!(
            (c.xi(&x, &pbc) - expected).abs() < TOL,
            "wrapped cyl ξ={}, expected {expected}",
            c.xi(&x, &pbc)
        );
        let raw = (8.0_f64 * 8.0 + 3.0 * 3.0).sqrt();
        assert!(
            (c.xi(&x, &PbcWrap::none()) - raw).abs() < TOL,
            "raw cyl ξ={}, expected {raw}",
            c.xi(&x, &PbcWrap::none())
        );
    }

    // ── construction validation ──────────────────────────────────────────────

    #[test]
    fn degenerate_normal_rejected() {
        assert_eq!(
            Coordinate::planar([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]).unwrap_err(),
            CoordinateError::ZeroNormal
        );
    }

    #[test]
    fn degenerate_axis_rejected() {
        assert_eq!(
            Coordinate::cylindrical([0.0; 3], [0.0, 0.0, 0.0], 0.1).unwrap_err(),
            CoordinateError::ZeroAxis
        );
    }

    #[test]
    fn constructors_normalize_directions() {
        // Planar normal of length 2 → ξ for a unit step equals the step, not 2×.
        let c = Coordinate::planar([2.0, 0.0, 0.0], [0.0; 3]).unwrap();
        assert!((c.xi(&[1.0, 0.0, 0.0], &PbcWrap::none()) - 1.0).abs() < TOL);
    }
}
