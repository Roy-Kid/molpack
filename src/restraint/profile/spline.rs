//! Tabulated target-distribution profile — monotone cubic spline + Boltzmann
//! inversion.
//!
//! This leaf module is the general C¹ distribution primitive underlying the
//! profile-distribution restraint family. It interpolates a user-supplied target
//! density ρ*(ξ) — a scalar reaction coordinate ξ (a length in Å) versus a number
//! density (count/volume) — tabulated on a strictly increasing grid, using a
//! monotone (Fritsch–Carlson clamped) cubic spline so no segment overshoots, and
//! performs numerical Boltzmann inversion to a soft energy `U(ξ)` and force
//! `−dU/dξ`:
//!
//! ```text
//! U(ξ)    = −kT·ln( max(s(ξ), ρ_min) / ρ₀ )
//! dU/dξ   = −kT·s'(ξ) / max(s(ξ), ρ_min)
//! ```
//!
//! where `s(ξ)` is the spline value, `s'(ξ)` its analytic derivative, and the
//! density floor `ρ_min` caps `U` so empty bins stay finite. The three closed-form
//! analytic families in [`distribution`](super::distribution) (Gaussian / erf /
//! exponential) are special cases: a tabulated sampling of any of them reproduces
//! the analytic `U` within interpolation tolerance.
//!
//! Units: `ξ` length (Å), `ρ*` number density (count/volume), `U`/`kT` energy.
//! `kT` is a FIXED energy-scale factor (NOT a simulated temperature) that sets the
//! bias steepness; `ρ₀` is an arbitrary reference density entering `U` only as an
//! additive constant, so it does not affect `dU/dξ` or the force.
//!
//! Construction-time concerns shared with [`distribution`](super::distribution)
//! are reused verbatim: the shell-volume [`ShellJacobian`] converts a count
//! histogram n(ξ) → volumetric density ρ*(ξ) = n(ξ)/(dV/dξ) BEFORE fitting (so
//! downstream evaluation always sees a density), [`InputKind`] flags whether the
//! supplied values are already a density, and [`DensityFloor`] supplies the floor
//! `ρ_min`, reference `ρ₀`, and energy cap.
//!
//! This module owns no coordinate geometry and no `Restraint` impl — it is a pure
//! numerical primitive consumed by a later sub-spec.

use molrs::types::F;

use super::distribution::{DensityFloor, InputKind, ShellJacobian};

/// Error returned by [`TabulatedProfile::new`] when the supplied grid or values
/// are not a valid tabulated density.
#[derive(Debug, Clone, PartialEq)]
pub enum SplineError {
    /// Fewer than two grid nodes were supplied (a spline needs ≥ 2).
    TooFewNodes(usize),
    /// The node count and value count differ.
    LengthMismatch {
        /// Number of grid nodes ξ_i.
        nodes: usize,
        /// Number of tabulated values.
        values: usize,
    },
    /// The grid ξ_i is not strictly increasing at the given index.
    NotStrictlyIncreasing(usize),
    /// A node ξ_i or value was non-finite (NaN / ∞) at the given index.
    NonFinite(usize),
    /// A tabulated value was negative at the given index (density must be ≥ 0).
    Negative(usize),
    /// A radial / cylindrical histogram grid included ξ ≤ 0, where the shell
    /// factor `dV/dξ` vanishes and the conversion would divide by zero.
    ShellSingularAtOrigin(usize),
}

impl std::fmt::Display for SplineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SplineError::TooFewNodes(n) => {
                write!(f, "tabulated profile needs >= 2 nodes, got {n}")
            }
            SplineError::LengthMismatch { nodes, values } => {
                write!(f, "node count {nodes} != value count {values}")
            }
            SplineError::NotStrictlyIncreasing(i) => {
                write!(f, "grid not strictly increasing at index {i}")
            }
            SplineError::NonFinite(i) => {
                write!(f, "non-finite node or value at index {i}")
            }
            SplineError::Negative(i) => {
                write!(f, "negative tabulated density at index {i}")
            }
            SplineError::ShellSingularAtOrigin(i) => write!(
                f,
                "radial/cylindrical histogram node at index {i} has xi <= 0 where dV/dxi vanishes"
            ),
        }
    }
}

impl std::error::Error for SplineError {}

/// A C¹ monotone (Fritsch–Carlson clamped) cubic spline of a tabulated target
/// density ρ*(ξ), with analytic derivative and Boltzmann inversion to `U(ξ)` /
/// `dU/dξ`.
///
/// Immutable, `Debug`, `Send + Sync`: it stores the sorted nodes, the volumetric
/// density at each node (after any histogram→density conversion), and the
/// per-segment Hermite slopes; no interior mutability. Outside `[ξ₀, ξ_last]`
/// evaluation uses flat extrapolation — `value` clamps to the end-node value and
/// `deriv` is zero — so value and derivative stay continuous at the boundary.
#[derive(Debug, Clone)]
pub struct TabulatedProfile {
    /// Strictly increasing grid nodes ξ_i (Å).
    nodes: Vec<F>,
    /// Volumetric density ρ*(ξ_i) at each node (after Jacobian conversion).
    rho: Vec<F>,
    /// Node derivatives m_i (the Fritsch–Carlson-limited Hermite slopes).
    slopes: Vec<F>,
    /// The density floor / energy cap shared with the analytic family.
    floor: DensityFloor,
}

impl TabulatedProfile {
    /// Build a profile from a grid `nodes`, raw tabulated `values`, an
    /// `input_kind`, the geometry's shell `jacobian`, and a density `floor`.
    ///
    /// Validates the grid (≥ 2 nodes, matched lengths, strictly increasing ξ, all
    /// nodes and values finite and non-negative). When `input_kind` is
    /// [`CountHistogram`](InputKind::CountHistogram) each value n(ξ_i) is converted
    /// to a volumetric density ρ*(ξ_i) = n(ξ_i)/(dV/dξ)(ξ_i) using `jacobian`
    /// BEFORE fitting, so downstream evaluation always sees a density. A
    /// radial/cylindrical histogram grid must not include ξ ≤ 0 (the shell factor
    /// `dV/dξ` vanishes there) — such a grid is rejected with
    /// [`ShellSingularAtOrigin`](SplineError::ShellSingularAtOrigin). The
    /// Fritsch–Carlson limiter then fixes the node slopes so the stored Hermite
    /// spline is monotone and overshoot-free.
    ///
    /// # Errors
    /// Returns a [`SplineError`] describing the first grid/value problem found.
    pub fn new(
        nodes: Vec<F>,
        values: Vec<F>,
        input_kind: InputKind,
        jacobian: ShellJacobian,
        floor: DensityFloor,
    ) -> Result<Self, SplineError> {
        if nodes.len() < 2 {
            return Err(SplineError::TooFewNodes(nodes.len()));
        }
        if nodes.len() != values.len() {
            return Err(SplineError::LengthMismatch {
                nodes: nodes.len(),
                values: values.len(),
            });
        }
        Self::validate(&nodes, &values)?;
        let rho = Self::to_density(&nodes, values, input_kind, jacobian)?;
        let slopes = fritsch_carlson_slopes(&nodes, &rho);
        Ok(Self {
            nodes,
            rho,
            slopes,
            floor,
        })
    }

    /// Validate finiteness, strict monotonicity of ξ, and non-negativity of the
    /// raw tabulated values.
    fn validate(nodes: &[F], values: &[F]) -> Result<(), SplineError> {
        for (i, (&x, &v)) in nodes.iter().zip(values.iter()).enumerate() {
            if !x.is_finite() || !v.is_finite() {
                return Err(SplineError::NonFinite(i));
            }
            if v < 0.0 {
                return Err(SplineError::Negative(i));
            }
            if i > 0 && x <= nodes[i - 1] {
                return Err(SplineError::NotStrictlyIncreasing(i));
            }
        }
        Ok(())
    }

    /// Convert raw `values` to a volumetric density: identity for a density input,
    /// or `n(ξ_i)/(dV/dξ)(ξ_i)` for a histogram input. Guards against the shell
    /// singularity at ξ ≤ 0 for radial/cylindrical geometries.
    fn to_density(
        nodes: &[F],
        values: Vec<F>,
        input_kind: InputKind,
        jacobian: ShellJacobian,
    ) -> Result<Vec<F>, SplineError> {
        if input_kind == InputKind::VolumetricDensity {
            return Ok(values);
        }
        let mut rho = values;
        for (i, (&x, v)) in nodes.iter().zip(rho.iter_mut()).enumerate() {
            let dv_dxi = shell_dv_dxi(jacobian, x);
            if dv_dxi <= 0.0 {
                return Err(SplineError::ShellSingularAtOrigin(i));
            }
            *v /= dv_dxi;
        }
        Ok(rho)
    }

    /// Locate the segment containing `xi`. Returns `Some(i)` for the segment
    /// `[ξ_i, ξ_{i+1}]`, or `None` when `xi` is outside `[ξ₀, ξ_last]`
    /// (flat-extrapolated region).
    #[inline]
    fn segment(&self, xi: F) -> Option<usize> {
        let last = self.nodes.len() - 1;
        if xi < self.nodes[0] || xi > self.nodes[last] {
            return None;
        }
        // Binary search for the rightmost node <= xi, clamped to a valid segment.
        let mut idx = self.nodes.partition_point(|&n| n <= xi).saturating_sub(1);
        if idx >= last {
            idx = last - 1;
        }
        Some(idx)
    }

    /// Spline value s(ξ): the monotone Hermite cubic on the bracketing segment,
    /// or the clamped end-node value under flat extrapolation outside the grid.
    pub fn value(&self, xi: F) -> F {
        match self.segment(xi) {
            None => {
                if xi < self.nodes[0] {
                    self.rho[0]
                } else {
                    self.rho[self.nodes.len() - 1]
                }
            }
            Some(i) => {
                let (h, t) = self.segment_param(i, xi);
                hermite_value(
                    self.rho[i],
                    self.rho[i + 1],
                    self.slopes[i],
                    self.slopes[i + 1],
                    h,
                    t,
                )
            }
        }
    }

    /// Analytic spline derivative s'(ξ); zero in the flat-extrapolated regions.
    pub fn deriv(&self, xi: F) -> F {
        match self.segment(xi) {
            None => 0.0,
            Some(i) => {
                let (h, t) = self.segment_param(i, xi);
                hermite_deriv(
                    self.rho[i],
                    self.rho[i + 1],
                    self.slopes[i],
                    self.slopes[i + 1],
                    h,
                    t,
                )
            }
        }
    }

    /// Segment width `h = ξ_{i+1} − ξ_i` and local parameter `t = (ξ − ξ_i)/h`.
    #[inline]
    fn segment_param(&self, i: usize, xi: F) -> (F, F) {
        let h = self.nodes[i + 1] - self.nodes[i];
        let t = (xi - self.nodes[i]) / h;
        (h, t)
    }

    /// Floored density `max(s(ξ), ρ_min)`.
    #[inline]
    fn floored(&self, xi: F) -> F {
        self.value(xi).max(self.floor.rho_min)
    }

    /// `U(ξ) = −kT·ln( max(s(ξ), ρ_min) / ρ₀ )` (energy). `ρ₀` enters only as an
    /// additive constant.
    pub fn u(&self, xi: F, kt: F) -> F {
        -kt * (self.floored(xi) / self.floor.rho0).ln()
    }

    /// `dU/dξ = −kT·s'(ξ)/max(s(ξ), ρ_min)` (force generator); zero where the
    /// floor is active (`s(ξ) ≤ ρ_min`) so a capped/empty bin exerts no force.
    pub fn du_dxi(&self, xi: F, kt: F) -> F {
        if self.value(xi) <= self.floor.rho_min {
            return 0.0;
        }
        -kt * self.deriv(xi) / self.floored(xi)
    }
}

/// The shell-volume Jacobian `dV/dξ` at `xi`: `1` planar (constant, no
/// centre-weighting), `4π ξ²` radial, `2π ξ L` cylindrical. Reuses the analytic
/// family's [`ShellJacobian`] via its `log_shell_correction` so the two paths
/// stay numerically identical.
#[inline]
fn shell_dv_dxi(jacobian: ShellJacobian, xi: F) -> F {
    match jacobian {
        ShellJacobian::Planar => 1.0,
        _ => jacobian.log_shell_correction(xi).exp(),
    }
}

/// Cubic Hermite value on a segment of width `h` at local parameter `t ∈ [0, 1]`,
/// given endpoint values `y0, y1` and endpoint slopes `m0, m1` (per unit ξ).
#[inline]
fn hermite_value(y0: F, y1: F, m0: F, m1: F, h: F, t: F) -> F {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
}

/// Analytic ξ-derivative of [`hermite_value`] (chain rule through `t = ·/h`).
#[inline]
fn hermite_deriv(y0: F, y1: F, m0: F, m1: F, h: F, t: F) -> F {
    let t2 = t * t;
    let dh00 = 6.0 * t2 - 6.0 * t;
    let dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    let dh01 = -6.0 * t2 + 6.0 * t;
    let dh11 = 3.0 * t2 - 2.0 * t;
    (dh00 * y0 + dh01 * y1) / h + dh10 * m0 + dh11 * m1
}

/// Fritsch–Carlson monotone node slopes for a Hermite spline of `(nodes, rho)`.
///
/// Computes secant slopes, seeds interior node slopes with a finite-difference
/// average (endpoints with the one-sided secant), then applies the Fritsch–Carlson
/// limiter so no segment overshoots: a flat secant forces a zero slope, and the
/// `(α, β)` pair is projected into the monotonicity circle of radius 3.
fn fritsch_carlson_slopes(nodes: &[F], rho: &[F]) -> Vec<F> {
    let n = nodes.len();
    // Secant slopes Δ_k for each segment [k, k+1].
    let deltas: Vec<F> = (0..n - 1)
        .map(|k| (rho[k + 1] - rho[k]) / (nodes[k + 1] - nodes[k]))
        .collect();

    // Initial node slopes m_i.
    let mut m = vec![0.0; n];
    m[0] = deltas[0];
    m[n - 1] = deltas[n - 2];
    for i in 1..n - 1 {
        // Zero slope at a local extremum or across a sign change avoids overshoot.
        if deltas[i - 1] * deltas[i] <= 0.0 {
            m[i] = 0.0;
        } else {
            m[i] = 0.5 * (deltas[i - 1] + deltas[i]);
        }
    }

    // Fritsch–Carlson limiter: project (α, β) into the monotone region.
    for k in 0..n - 1 {
        if deltas[k] == 0.0 {
            // Flat segment ⇒ both bounding slopes must be zero (no overshoot).
            m[k] = 0.0;
            m[k + 1] = 0.0;
            continue;
        }
        let alpha = m[k] / deltas[k];
        let beta = m[k + 1] / deltas[k];
        let s = alpha * alpha + beta * beta;
        if s > 9.0 {
            let tau = 3.0 / s.sqrt();
            m[k] = tau * alpha * deltas[k];
            m[k + 1] = tau * beta * deltas[k];
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const KT: F = 0.6;

    fn floor(rho_min: F) -> DensityFloor {
        DensityFloor::new(rho_min, 1.0, KT)
    }

    /// Build a density-input profile (planar geometry, no Jacobian) from
    /// matched node/value vectors with a tiny floor.
    fn density_profile(nodes: Vec<F>, values: Vec<F>) -> TabulatedProfile {
        TabulatedProfile::new(
            nodes,
            values,
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .expect("valid grid")
    }

    fn assert_close(got: F, want: F, tol: F, label: &str) {
        assert!(
            (got - want).abs() <= tol,
            "{label}: got={got} want={want} diff={}",
            (got - want).abs()
        );
    }

    #[test]
    fn profile_is_send_sync_debug() {
        fn assert_send_sync<T: Send + Sync + std::fmt::Debug>() {}
        assert_send_sync::<TabulatedProfile>();
        let p = density_profile(vec![0.0, 1.0], vec![1.0, 2.0]);
        let _ = format!("{p:?}");
    }

    // ── grid validation ───────────────────────────────────────────────────────

    #[test]
    fn rejects_too_few_nodes() {
        let err = TabulatedProfile::new(
            vec![1.0],
            vec![2.0],
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(err, SplineError::TooFewNodes(1));
    }

    #[test]
    fn rejects_length_mismatch() {
        let err = TabulatedProfile::new(
            vec![1.0, 2.0, 3.0],
            vec![1.0, 2.0],
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(
            err,
            SplineError::LengthMismatch {
                nodes: 3,
                values: 2
            }
        );
    }

    #[test]
    fn rejects_non_increasing_grid() {
        let err = TabulatedProfile::new(
            vec![1.0, 1.0, 2.0],
            vec![1.0, 2.0, 3.0],
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(err, SplineError::NotStrictlyIncreasing(1));
    }

    #[test]
    fn rejects_non_finite_value() {
        let err = TabulatedProfile::new(
            vec![1.0, 2.0, 3.0],
            vec![1.0, F::NAN, 3.0],
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(err, SplineError::NonFinite(1));
    }

    #[test]
    fn rejects_negative_value() {
        let err = TabulatedProfile::new(
            vec![1.0, 2.0, 3.0],
            vec![1.0, -0.5, 3.0],
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(err, SplineError::Negative(1));
    }

    #[test]
    fn rejects_radial_histogram_at_origin() {
        let err = TabulatedProfile::new(
            vec![0.0, 1.0, 2.0],
            vec![1.0, 1.0, 1.0],
            InputKind::CountHistogram,
            ShellJacobian::Radial,
            floor(1e-30),
        )
        .unwrap_err();
        assert_eq!(err, SplineError::ShellSingularAtOrigin(0));
    }

    // ── ac-001: exact node interpolation + C¹ continuity ───────────────────────

    #[test]
    fn passes_through_every_node() {
        let nodes = vec![0.0, 1.0, 2.5, 4.0, 5.0, 7.0];
        let values = vec![0.2, 1.3, 0.7, 2.1, 1.9, 0.4];
        let p = density_profile(nodes.clone(), values.clone());
        for (xi, want) in nodes.iter().zip(values.iter()) {
            assert_close(p.value(*xi), *want, 1e-12, &format!("node @ {xi}"));
        }
    }

    #[test]
    fn c1_continuity_across_interior_nodes() {
        let nodes = vec![0.0, 1.0, 2.5, 4.0, 5.0, 7.0];
        let values = vec![0.2, 1.3, 0.7, 2.1, 1.9, 0.4];
        let p = density_profile(nodes.clone(), values.clone());
        // Approach each interior node from the left and right segments; value and
        // derivative must match (C¹).
        let eps = 1e-7;
        for &xi in &nodes[1..nodes.len() - 1] {
            let vl = p.value(xi - eps);
            let vr = p.value(xi + eps);
            assert_close(vl, vr, 1e-5, &format!("value continuity @ {xi}"));
            let dl = p.deriv(xi - eps);
            let dr = p.deriv(xi + eps);
            assert_close(dl, dr, 1e-4, &format!("deriv continuity @ {xi}"));
        }
    }

    // ── ac-002: analytic deriv vs central FD of value ──────────────────────────

    #[test]
    fn analytic_deriv_matches_central_fd() {
        let nodes = vec![0.0, 1.0, 2.5, 4.0, 5.0, 7.0];
        let values = vec![0.2, 1.3, 0.7, 2.1, 1.9, 0.4];
        let p = density_profile(nodes.clone(), values.clone());
        let h = 1e-6;
        // Sample the interior, avoiding nodes (where the FD straddles a curvature
        // jump) and the flat extrapolated regions.
        let mut xi = nodes[0] + 0.05;
        while xi < nodes[nodes.len() - 1] - 0.05 {
            // Skip a small neighborhood of each node.
            let near_node = nodes.iter().any(|n| (xi - n).abs() < 2.0 * h);
            if !near_node {
                let fd = (p.value(xi + h) - p.value(xi - h)) / (2.0 * h);
                assert_close(p.deriv(xi), fd, 1e-4, &format!("deriv vs fd @ {xi}"));
            }
            xi += 0.1;
        }
    }

    // ── flat extrapolation outside the grid ────────────────────────────────────

    #[test]
    fn flat_extrapolation_outside_grid() {
        let nodes = vec![1.0, 2.0, 3.0];
        let values = vec![0.5, 1.5, 0.8];
        let p = density_profile(nodes.clone(), values.clone());
        // Below first node: end-node value, zero derivative.
        assert_close(p.value(0.0), 0.5, 1e-12, "below value");
        assert_close(p.deriv(0.0), 0.0, 1e-12, "below deriv");
        // Above last node.
        assert_close(p.value(5.0), 0.8, 1e-12, "above value");
        assert_close(p.deriv(5.0), 0.0, 1e-12, "above deriv");
        // Continuity at the boundary nodes (one-sided approach from inside).
        let eps = 1e-7;
        assert_close(p.value(1.0 + eps), 0.5, 1e-5, "left-boundary continuity");
        assert_close(p.value(3.0 - eps), 0.8, 1e-5, "right-boundary continuity");
    }

    // ── ac-004(a): monotone input is overshoot-free ────────────────────────────

    #[test]
    fn monotone_decreasing_input_has_no_overshoot() {
        // Steep monotone-decreasing data: a naive natural cubic overshoots.
        let nodes = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let values = vec![10.0, 9.8, 9.5, 1.0, 0.5, 0.2];
        let p = density_profile(nodes.clone(), values.clone());
        // Sample densely; every value must lie within the bracket of its segment's
        // endpoints (no overshoot above the upper node or below the lower node).
        for i in 0..nodes.len() - 1 {
            let (lo_x, hi_x) = (nodes[i], nodes[i + 1]);
            let (a, b) = (values[i], values[i + 1]);
            let (bracket_lo, bracket_hi) = (a.min(b), a.max(b));
            let mut t = lo_x;
            while t <= hi_x {
                let v = p.value(t);
                assert!(
                    v >= bracket_lo - 1e-9 && v <= bracket_hi + 1e-9,
                    "overshoot in segment [{lo_x},{hi_x}] @ {t}: v={v} bracket=[{bracket_lo},{bracket_hi}]"
                );
                t += 0.02;
            }
        }
    }

    // ── ac-004(b): zero bin yields finite floored u / du ───────────────────────

    #[test]
    fn zero_bin_is_finite_via_floor() {
        let rho_min = 1e-3;
        let nodes = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![1.0, 0.0, 0.0, 1.0]; // empty middle bins
        let p = TabulatedProfile::new(
            nodes,
            values,
            InputKind::VolumetricDensity,
            ShellJacobian::Planar,
            floor(rho_min),
        )
        .expect("valid");
        let u_max = -KT * (rho_min / 1.0).ln();
        for &xi in &[1.0, 1.5, 2.0] {
            let u = p.u(xi, KT);
            let du = p.du_dxi(xi, KT);
            assert!(u.is_finite(), "u not finite @ {xi}: {u}");
            assert!(du.is_finite(), "du not finite @ {xi}: {du}");
            // At the floored minimum, U is capped and force is zero.
            assert!(u <= u_max + 1e-9, "u exceeds cap @ {xi}: {u} > {u_max}");
        }
        // Exactly at the empty node the density is 0 -> capped, no force.
        assert_eq!(p.du_dxi(1.5, KT), 0.0, "floored bin must exert no force");
    }

    // ── ac-003: tabulated Gaussian reproduces harmonic U ───────────────────────

    #[test]
    fn tabulated_gaussian_reproduces_harmonic_u() {
        let (mu, sigma) = (6.0, 1.5);
        // Sample rho* ∝ exp(−(ξ−μ)²/2σ²) onto a fine grid spanning ±3σ.
        let n = 61;
        let (lo, hi) = (mu - 3.0 * sigma, mu + 3.0 * sigma);
        let nodes: Vec<F> = (0..n)
            .map(|i| lo + (hi - lo) * (i as F) / ((n - 1) as F))
            .collect();
        let gauss = |xi: F| (-0.5 * ((xi - mu) / sigma).powi(2)).exp();
        let values: Vec<F> = nodes.iter().map(|&x| gauss(x)).collect();
        let p = density_profile(nodes.clone(), values);
        // Compare U DIFFERENCES against the harmonic well (the additive ρ₀ const
        // cancels). U(ξ) − U(μ) = (kT/2σ²)(ξ−μ)².
        let u_mu = p.u(mu, KT);
        let k = KT / (2.0 * sigma * sigma);
        for &xi in &[3.5, 4.5, 5.5, 6.0, 6.5, 7.5, 8.5] {
            let want = k * (xi - mu).powi(2);
            let got = p.u(xi, KT) - u_mu;
            assert_close(got, want, 5e-3, &format!("harmonic U @ {xi}"));
        }
    }

    // ── ac-004(c): radial histogram needs the ∝ξ² Jacobian ─────────────────────

    #[test]
    fn radial_histogram_recovers_uniform_density_only_with_jacobian() {
        // A uniform volumetric density rho* = C means the count histogram is
        // n(ξ) = C·4πξ². Feed that n(ξ) as a CountHistogram+Radial profile; after
        // the Jacobian conversion the stored density must be ~C (flat).
        let c = 0.5;
        let nodes = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let counts: Vec<F> = nodes.iter().map(|&x| c * 4.0 * PI * x * x).collect();
        let with_jac = TabulatedProfile::new(
            nodes.clone(),
            counts.clone(),
            InputKind::CountHistogram,
            ShellJacobian::Radial,
            floor(1e-30),
        )
        .expect("valid");
        // The recovered density is flat at C across the grid.
        for &xi in &nodes {
            assert_close(with_jac.value(xi), c, 1e-9, &format!("recovered C @ {xi}"));
        }
        // Without the Jacobian (treat n as a density), the stored value is the raw
        // count C·4πξ² — strongly ξ-dependent, NOT flat.
        let without_jac = density_profile(nodes.clone(), counts.clone());
        let v1 = without_jac.value(1.0);
        let v5 = without_jac.value(5.0);
        assert!(
            (v5 / v1 - 25.0).abs() < 1e-6,
            "without Jacobian density should scale as xi^2 (ratio ~25), got {}",
            v5 / v1
        );
    }
}
