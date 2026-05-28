//! Exact port of `polartocart.f90` from Packmol.
//!
//! Euler angle convention (eulerrmat). Packmol spells these `beta`/`gama`/
//! `teta`; renamed here for clarity (`gama`/`teta` were Fortran misspellings):
//!   euler_beta  = rotation about y-axis  (Packmol: `beta`)
//!   euler_gamma = rotation about z-axis  (Packmol: `gama`)
//!   euler_theta = rotation about x-axis  (Packmol: `teta`)

use molrs::types::F;
/// Compute rotation matrix columns from Euler angles.
/// Port of Fortran `eulerrmat`.
///
/// Returns (v1, v2, v3) — the three columns of the rotation matrix.
#[inline(always)]
pub fn eulerrmat(euler_beta: F, euler_gamma: F, euler_theta: F) -> ([F; 3], [F; 3], [F; 3]) {
    let cb = euler_beta.cos();
    let sb = euler_beta.sin();
    let cg = euler_gamma.cos();
    let sg = euler_gamma.sin();
    let ct = euler_theta.cos();
    let st = euler_theta.sin();

    let v1 = [-sb * sg * ct + cb * cg, -sb * cg * ct - cb * sg, sb * st];
    let v2 = [cb * sg * ct + sb * cg, cb * cg * ct - sb * sg, -cb * st];
    let v3 = [sg * st, cg * st, ct];

    (v1, v2, v3)
}

/// Compute Cartesian coordinates from center-of-mass, reference coordinates, and rotation matrix.
/// Port of Fortran `compcart`.
#[inline(always)]
pub fn compcart(xcm: &[F; 3], xref: &[F; 3], v1: &[F; 3], v2: &[F; 3], v3: &[F; 3]) -> [F; 3] {
    [
        xcm[0] + xref[0] * v1[0] + xref[1] * v2[0] + xref[2] * v3[0],
        xcm[1] + xref[0] * v1[1] + xref[1] * v2[1] + xref[2] * v3[1],
        xcm[2] + xref[0] * v1[2] + xref[1] * v2[2] + xref[2] * v3[2],
    ]
}

/// Compute rotation matrix for "fixed" molecules using the "human" convention.
/// Port of Fortran `eulerfixed`.
///
/// In this convention:
///   euler_beta  = counterclockwise rotation around x-axis
///   euler_gamma  = counterclockwise rotation around y-axis
///   euler_theta  = counterclockwise rotation around z-axis
#[inline(always)]
pub fn eulerfixed(euler_beta: F, euler_gamma: F, euler_theta: F) -> ([F; 3], [F; 3], [F; 3]) {
    let c1 = euler_beta.cos();
    let s1 = euler_beta.sin();
    let c2 = euler_gamma.cos();
    let s2 = euler_gamma.sin();
    let c3 = euler_theta.cos();
    let s3 = euler_theta.sin();

    let v1 = [c2 * c3, c1 * s3 + c3 * s1 * s2, s1 * s3 - c1 * c3 * s2];
    let v2 = [-c2 * s3, c1 * c3 - s1 * s2 * s3, c1 * s2 * s3 + c3 * s1];
    let v3 = [s2, -c2 * s1, c1 * c2];

    (v1, v2, v3)
}

/// All 9 partial derivatives of rotation matrix columns w.r.t. euler_beta/euler_gamma/euler_theta.
/// Exact port of `computeg.f90` lines 169-204.
///
/// Returns (dv1beta, dv1gama, dv1teta, dv2beta, dv2gama, dv2teta, dv3beta, dv3gama, dv3teta)
#[allow(clippy::type_complexity)]
pub fn eulerrmat_derivatives(
    euler_beta: F,
    euler_gamma: F,
    euler_theta: F,
) -> (
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
    [F; 3],
) {
    let cb = euler_beta.cos();
    let sb = euler_beta.sin();
    let cg = euler_gamma.cos();
    let sg = euler_gamma.sin();
    let ct = euler_theta.cos();
    let st = euler_theta.sin();

    let dv1beta = [-cb * sg * ct - sb * cg, -cb * cg * ct + sb * sg, cb * st];
    let dv2beta = [-sb * sg * ct + cb * cg, -sb * cg * ct - cb * sg, sb * st];
    let dv3beta = [0.0, 0.0, 0.0];

    let dv1gama = [-sb * cg * ct - cb * sg, sb * sg * ct - cb * cg, 0.0];
    let dv2gama = [cb * cg * ct - sb * sg, -sg * cb * ct - cg * sb, 0.0];
    let dv3gama = [cg * st, -sg * st, 0.0];

    let dv1teta = [sb * sg * st, sb * cg * st, sb * ct];
    let dv2teta = [-cb * sg * st, -cb * cg * st, -cb * ct];
    let dv3teta = [sg * ct, cg * ct, -st];

    (
        dv1beta, dv1gama, dv1teta, dv2beta, dv2gama, dv2teta, dv3beta, dv3gama, dv3teta,
    )
}
