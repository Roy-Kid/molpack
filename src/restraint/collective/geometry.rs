//! Reaction-coordinate geometries for the distribution-matching family, plus the
//! thin composition that turns a geometry + a target quantile function into the
//! `f` / `fg` a [`Restraint`](super::Restraint) needs.
//!
//! A geometry maps a group's Cartesian coordinates to a scalar ξ and scatters
//! `∂L/∂ξ` back onto them via `∇ξ`:
//!
//! - **plane** — ξ = x·n̂ − offset (signed distance to a plane); `∇ξ = n̂`.
//! - **point** — ξ = ‖x − c‖ (distance to a centre); `∇ξ = (x − c)/‖x − c‖`,
//!   guarded at the centre (where it is undefined) so no spurious force results.
//!
//! Each geometry then exposes `*_match_f` / `*_match_fg`: compute ξ, run the
//! shared Wasserstein [`engine`](super::engine) against the supplied `quantile`,
//! and (for `fg`) scatter the coupled gradient. A concrete `<Distribution><Geometry>`
//! type is therefore a one-line delegation per method — adding a geometry is a
//! new `*_xi`/`*_scatter` pair plus its two `*_match_*` wrappers, reused by every
//! distribution.

use molrs::types::F;

use super::engine::{wasserstein_grad, wasserstein_value};

// --------------------------------------------------------------------- plane

/// ξᵢ = xᵢ·n̂ − offset for a plane with **unit** normal `normal`.
fn plane_xi(coords: &[[F; 3]], normal: &[F; 3], offset: F) -> Vec<F> {
    coords
        .iter()
        .map(|x| x[0] * normal[0] + x[1] * normal[1] + x[2] * normal[2] - offset)
        .collect()
}

/// Accumulate `∂L/∂ξᵢ · n̂` into `grads[i]` (∇ξ is the constant plane normal).
fn plane_scatter(dxi: &[F], normal: &[F; 3], grads: &mut [[F; 3]]) {
    for (g, &s) in grads.iter_mut().zip(dxi.iter()) {
        g[0] += s * normal[0];
        g[1] += s * normal[1];
        g[2] += s * normal[2];
    }
}

/// Value of the distribution-match penalty for the plane coordinate.
pub(super) fn plane_match_f(
    coords: &[[F; 3]],
    normal: &[F; 3],
    offset: F,
    strength: F,
    quantile: impl Fn(F) -> F,
) -> F {
    wasserstein_value(&plane_xi(coords, normal, offset), strength, quantile)
}

/// Value + gradient of the distribution-match penalty for the plane coordinate;
/// the coupled gradient is scattered into `grads`.
pub(super) fn plane_match_fg(
    coords: &[[F; 3]],
    normal: &[F; 3],
    offset: F,
    strength: F,
    quantile: impl Fn(F) -> F,
    grads: &mut [[F; 3]],
) -> F {
    let xi = plane_xi(coords, normal, offset);
    let mut dxi = vec![0.0 as F; xi.len()];
    let e = wasserstein_grad(&xi, strength, quantile, &mut dxi);
    plane_scatter(&dxi, normal, grads);
    e
}

// --------------------------------------------------------------------- point

/// ξᵢ = ‖xᵢ − c‖ (distance from the centre `center`).
fn point_xi(coords: &[[F; 3]], center: &[F; 3]) -> Vec<F> {
    coords
        .iter()
        .map(|x| {
            let d = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .collect()
}

/// Below this radius the outward unit vector r̂ is treated as undefined; the
/// atom contributes no radial force (it sits on the centre, a measure-zero case).
const R_GUARD: F = 1e-9;

/// Accumulate `∂L/∂ξᵢ · r̂ᵢ` into `grads[i]`, where `r̂ᵢ = (xᵢ − c)/‖xᵢ − c‖`.
fn point_scatter(dxi: &[F], coords: &[[F; 3]], center: &[F; 3], grads: &mut [[F; 3]]) {
    for (i, g) in grads.iter_mut().enumerate() {
        let d = [
            coords[i][0] - center[0],
            coords[i][1] - center[1],
            coords[i][2] - center[2],
        ];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        if r > R_GUARD {
            let s = dxi[i] / r;
            g[0] += s * d[0];
            g[1] += s * d[1];
            g[2] += s * d[2];
        }
    }
}

/// Value of the distribution-match penalty for the point (radial) coordinate.
pub(super) fn point_match_f(
    coords: &[[F; 3]],
    center: &[F; 3],
    strength: F,
    quantile: impl Fn(F) -> F,
) -> F {
    wasserstein_value(&point_xi(coords, center), strength, quantile)
}

/// Value + gradient of the distribution-match penalty for the point coordinate;
/// the coupled gradient is scattered into `grads`.
pub(super) fn point_match_fg(
    coords: &[[F; 3]],
    center: &[F; 3],
    strength: F,
    quantile: impl Fn(F) -> F,
    grads: &mut [[F; 3]],
) -> F {
    let xi = point_xi(coords, center);
    let mut dxi = vec![0.0 as F; xi.len()];
    let e = wasserstein_grad(&xi, strength, quantile, &mut dxi);
    point_scatter(&dxi, coords, center, grads);
    e
}
