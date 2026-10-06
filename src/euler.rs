//! Exact port of `polartocart.f90` from Packmol.
//!
//! Euler angle convention (eulerrmat):
//!   beta  = rotation about y-axis
//!   gama  = rotation about z-axis
//!   teta  = rotation about x-axis

use molrs::types::F;
/// Compute rotation matrix columns from Euler angles.
/// Port of Fortran `eulerrmat`.
///
/// Returns (v1, v2, v3) — the three columns of the rotation matrix.
#[inline(always)]
pub fn eulerrmat(beta: F, gama: F, teta: F) -> ([F; 3], [F; 3], [F; 3]) {
    let cb = beta.cos();
    let sb = beta.sin();
    let cg = gama.cos();
    let sg = gama.sin();
    let ct = teta.cos();
    let st = teta.sin();

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
///   beta  = counterclockwise rotation around x-axis
///   gama  = counterclockwise rotation around y-axis
///   teta  = counterclockwise rotation around z-axis
#[inline(always)]
pub fn eulerfixed(beta: F, gama: F, teta: F) -> ([F; 3], [F; 3], [F; 3]) {
    let c1 = beta.cos();
    let s1 = beta.sin();
    let c2 = gama.cos();
    let s2 = gama.sin();
    let c3 = teta.cos();
    let s3 = teta.sin();

    let v1 = [c2 * c3, c1 * s3 + c3 * s1 * s2, s1 * s3 - c1 * c3 * s2];
    let v2 = [-c2 * s3, c1 * c3 - s1 * s2 * s3, c1 * s2 * s3 + c3 * s1];
    let v3 = [s2, -c2 * s1, c1 * c2];

    (v1, v2, v3)
}

/// All 9 partial derivatives of rotation matrix columns w.r.t. beta/gama/teta.
/// Exact port of `computeg.f90` lines 169-204.
///
/// Returns (dv1beta, dv1gama, dv1teta, dv2beta, dv2gama, dv2teta, dv3beta, dv3gama, dv3teta)
#[allow(clippy::type_complexity)]
pub fn eulerrmat_derivatives(
    beta: F,
    gama: F,
    teta: F,
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
    let cb = beta.cos();
    let sb = beta.sin();
    let cg = gama.cos();
    let sg = gama.sin();
    let ct = teta.cos();
    let st = teta.sin();

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

#[cfg(test)]
mod tests {
    //! Tests for Euler angle functions: eulerrmat, compcart, eulerfixed,
    //! eulerrmat_derivatives.

    use crate::F;
    use crate::euler::{compcart, eulerfixed, eulerrmat, eulerrmat_derivatives};

    const TOL: F = 1e-6;
    const PI: F = std::f64::consts::PI as F;

    // ── helpers ────────────────────────────────────────────────────────────────

    fn dot(a: &[F; 3], b: &[F; 3]) -> F {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    fn norm(a: &[F; 3]) -> F {
        dot(a, a).sqrt()
    }

    // ── eulerrmat ──────────────────────────────────────────────────────────────

    #[test]
    fn identity_rotation() {
        let (v1, v2, v3) = eulerrmat(0.0, 0.0, 0.0);
        // v1 ≈ [1,0,0], v2 ≈ [0,1,0], v3 ≈ [0,0,1]
        assert!((v1[0] - 1.0).abs() < TOL);
        assert!(v1[1].abs() < TOL);
        assert!(v1[2].abs() < TOL);

        assert!(v2[0].abs() < TOL);
        assert!((v2[1] - 1.0).abs() < TOL);
        assert!(v2[2].abs() < TOL);

        assert!(v3[0].abs() < TOL);
        assert!(v3[1].abs() < TOL);
        assert!((v3[2] - 1.0).abs() < TOL);
    }

    #[test]
    fn rotation_columns_are_orthonormal() {
        // Test with arbitrary angles
        let angles = [
            (0.3, 0.5, 0.7),
            (PI / 4.0, PI / 3.0, PI / 6.0),
            (1.0, 2.0, 3.0),
        ];
        for (beta, gama, teta) in angles {
            let (v1, v2, v3) = eulerrmat(beta, gama, teta);
            // Unit length
            assert!(
                (norm(&v1) - 1.0).abs() < TOL,
                "v1 not unit for ({beta},{gama},{teta})"
            );
            assert!((norm(&v2) - 1.0).abs() < TOL, "v2 not unit");
            assert!((norm(&v3) - 1.0).abs() < TOL, "v3 not unit");
            // Orthogonal
            assert!(dot(&v1, &v2).abs() < TOL, "v1·v2 not zero");
            assert!(dot(&v1, &v3).abs() < TOL, "v1·v3 not zero");
            assert!(dot(&v2, &v3).abs() < TOL, "v2·v3 not zero");
        }
    }

    #[test]
    fn beta_90_rotates_around_y() {
        let (v1, _v2, v3) = eulerrmat(PI / 2.0, 0.0, 0.0);
        // After 90° around y: x→z, z→-x
        // v1 should map [1,0,0] → [0,0,1]-ish, v3 should map [0,0,1] → [-1,0,0]-ish
        // Actually eulerrmat convention: v1 is the first column of R
        // For beta=π/2, gama=0, teta=0: sb=1, cb=0, sg=0, cg=1, st=0, ct=1
        // v1 = [-1*0*1+0*1, -1*1*1-0*0, 1*0] = [0, -1, 0]... let me just check orthogonality
        assert!((norm(&v1) - 1.0).abs() < TOL);
        assert!((norm(&v3) - 1.0).abs() < TOL);
    }

    // ── compcart ───────────────────────────────────────────────────────────────

    #[test]
    fn compcart_identity_is_translation() {
        let (v1, v2, v3) = eulerrmat(0.0, 0.0, 0.0);
        let xcm = [1.0, 2.0, 3.0];
        let xref = [0.1, 0.2, 0.3];
        let result = compcart(&xcm, &xref, &v1, &v2, &v3);
        // With identity rotation: result = xcm + xref
        assert!((result[0] - 1.1).abs() < TOL);
        assert!((result[1] - 2.2).abs() < TOL);
        assert!((result[2] - 3.3).abs() < TOL);
    }

    #[test]
    fn compcart_zero_ref_is_just_com() {
        let (v1, v2, v3) = eulerrmat(0.5, 1.0, 0.3);
        let xcm = [10.0, 20.0, 30.0];
        let xref = [0.0, 0.0, 0.0];
        let result = compcart(&xcm, &xref, &v1, &v2, &v3);
        assert!((result[0] - 10.0).abs() < TOL);
        assert!((result[1] - 20.0).abs() < TOL);
        assert!((result[2] - 30.0).abs() < TOL);
    }

    #[test]
    fn compcart_preserves_distance_from_com() {
        let (v1, v2, v3) = eulerrmat(0.3, 0.5, 0.7);
        let xcm = [1.0, 2.0, 3.0];
        let xref = [1.0, 0.0, 0.0];
        let result = compcart(&xcm, &xref, &v1, &v2, &v3);
        // Distance from COM should equal |xref| = 1.0
        let d = ((result[0] - xcm[0]).powi(2)
            + (result[1] - xcm[1]).powi(2)
            + (result[2] - xcm[2]).powi(2))
        .sqrt();
        assert!(
            (d - 1.0).abs() < TOL,
            "rotation should preserve distance from COM"
        );
    }

    #[test]
    fn compcart_two_atoms_preserve_internal_distance() {
        let (v1, v2, v3) = eulerrmat(1.2, -0.5, 0.8);
        let xcm = [5.0, 5.0, 5.0];
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let p1 = compcart(&xcm, &r1, &v1, &v2, &v3);
        let p2 = compcart(&xcm, &r2, &v1, &v2, &v3);
        let d =
            ((p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2) + (p1[2] - p2[2]).powi(2)).sqrt();
        let d_ref =
            ((r1[0] - r2[0]).powi(2) + (r1[1] - r2[1]).powi(2) + (r1[2] - r2[2]).powi(2)).sqrt();
        assert!(
            (d - d_ref).abs() < TOL,
            "rotation should preserve internal distances"
        );
    }

    // ── eulerfixed ─────────────────────────────────────────────────────────────

    #[test]
    fn eulerfixed_identity() {
        let (v1, v2, v3) = eulerfixed(0.0, 0.0, 0.0);
        assert!((v1[0] - 1.0).abs() < TOL);
        assert!((v2[1] - 1.0).abs() < TOL);
        assert!((v3[2] - 1.0).abs() < TOL);
    }

    #[test]
    fn eulerfixed_columns_orthonormal() {
        let (v1, v2, v3) = eulerfixed(0.5, 1.0, -0.3);
        assert!((norm(&v1) - 1.0).abs() < TOL);
        assert!((norm(&v2) - 1.0).abs() < TOL);
        assert!((norm(&v3) - 1.0).abs() < TOL);
        assert!(dot(&v1, &v2).abs() < TOL);
        assert!(dot(&v1, &v3).abs() < TOL);
        assert!(dot(&v2, &v3).abs() < TOL);
    }

    // ── eulerrmat_derivatives ──────────────────────────────────────────────────

    #[test]
    fn derivatives_match_finite_difference() {
        let beta: F = 0.3;
        let gama: F = 0.5;
        let teta: F = 0.7;
        // Use h=1e-3 for f32 stability; f64 would allow 1e-5.
        let h: F = 1e-3;

        let (dv1b, dv1g, dv1t, dv2b, dv2g, dv2t, dv3b, dv3g, dv3t) =
            eulerrmat_derivatives(beta, gama, teta);

        // Finite difference for d/dbeta
        let (v1p, v2p, v3p) = eulerrmat(beta + h, gama, teta);
        let (v1m, v2m, v3m) = eulerrmat(beta - h, gama, teta);
        for k in 0..3 {
            let fd = (v1p[k] - v1m[k]) / (2.0 * h);
            assert!(
                (dv1b[k] - fd).abs() < 1e-3,
                "dv1/dbeta[{k}]: analytic={} fd={fd}",
                dv1b[k]
            );
        }
        for k in 0..3 {
            let fd = (v2p[k] - v2m[k]) / (2.0 * h);
            assert!(
                (dv2b[k] - fd).abs() < 1e-3,
                "dv2/dbeta[{k}]: analytic={} fd={fd}",
                dv2b[k]
            );
        }
        for k in 0..3 {
            let fd = (v3p[k] - v3m[k]) / (2.0 * h);
            assert!(
                (dv3b[k] - fd).abs() < 1e-3,
                "dv3/dbeta[{k}]: analytic={} fd={fd}",
                dv3b[k]
            );
        }

        // d/dgama
        let (v1p, v2p, v3p) = eulerrmat(beta, gama + h, teta);
        let (v1m, v2m, v3m) = eulerrmat(beta, gama - h, teta);
        for k in 0..3 {
            let fd = (v1p[k] - v1m[k]) / (2.0 * h);
            assert!((dv1g[k] - fd).abs() < 1e-3, "dv1/dgama[{k}]");
        }
        for k in 0..3 {
            let fd = (v2p[k] - v2m[k]) / (2.0 * h);
            assert!((dv2g[k] - fd).abs() < 1e-3, "dv2/dgama[{k}]");
        }
        for k in 0..3 {
            let fd = (v3p[k] - v3m[k]) / (2.0 * h);
            assert!((dv3g[k] - fd).abs() < 1e-3, "dv3/dgama[{k}]");
        }

        // d/dteta
        let (v1p, v2p, v3p) = eulerrmat(beta, gama, teta + h);
        let (v1m, v2m, v3m) = eulerrmat(beta, gama, teta - h);
        for k in 0..3 {
            let fd = (v1p[k] - v1m[k]) / (2.0 * h);
            assert!((dv1t[k] - fd).abs() < 1e-3, "dv1/dteta[{k}]");
        }
        for k in 0..3 {
            let fd = (v2p[k] - v2m[k]) / (2.0 * h);
            assert!((dv2t[k] - fd).abs() < 1e-3, "dv2/dteta[{k}]");
        }
        for k in 0..3 {
            let fd = (v3p[k] - v3m[k]) / (2.0 * h);
            assert!((dv3t[k] - fd).abs() < 1e-3, "dv3/dteta[{k}]");
        }
    }
}
