//! Unit + property tests for the composed [`ProfileRestraint`] (sub-spec -05).
//!
//! Covers all four acceptance criteria:
//! - ac-001: planar-Gaussian `f` closed form; `fg` additive into a seeded `g`.
//! - ac-002: the full 4×5 coordinate×distribution gradient matrix vs central FD.
//! - ac-003: guard + floor keep `f`/`g` finite at a singular site / empty bin.
//! - ac-004: `Send + Sync + Debug`, no interior mutability, linear in `scale`.

use std::sync::Arc;

use molrs::types::F;

use crate::region::InsideBoxRegion;
use crate::restraint::Restraint;
use crate::restraint::profile::{
    Coordinate, DensityFloor, Distribution, InputKind, PbcWrap, ProfilePenalty, ProfileRestraint,
    ProfileTarget, ShellJacobian, TabulatedProfile,
};

const KT: F = 0.6;
const TOL: F = 1e-6;

/// A planar [`ProfilePenalty`] from one of the analytic distribution shapes,
/// with the safe defaults: no shell Jacobian (planar), volumetric density, and a
/// tiny floor that never trips for the in-range samples used here.
fn planar_penalty(dist: Distribution) -> ProfilePenalty {
    ProfilePenalty {
        dist,
        jacobian: ShellJacobian::Planar,
        input_kind: InputKind::VolumetricDensity,
        floor: DensityFloor::new(1e-30, 1.0, KT),
    }
}

/// A tabulated profile sampling a smooth Gaussian density on a fine grid — a
/// stand-in tabulated target whose derivative is smooth in the sampled interior.
fn tabulated_gaussian(mu: F, sigma: F) -> TabulatedProfile {
    let n = 81;
    let (lo, hi) = (mu - 4.0 * sigma, mu + 4.0 * sigma);
    let nodes: Vec<F> = (0..n)
        .map(|i| lo + (hi - lo) * (i as F) / ((n - 1) as F))
        .collect();
    let values: Vec<F> = nodes
        .iter()
        .map(|&x| (-0.5 * ((x - mu) / sigma).powi(2)).exp())
        .collect();
    TabulatedProfile::new(
        nodes,
        values,
        InputKind::VolumetricDensity,
        ShellJacobian::Planar,
        DensityFloor::new(1e-30, 1.0, KT),
    )
    .expect("valid tabulated grid")
}

// ── ac-001: planar-Gaussian closed form; fg additive ────────────────────────

#[test]
fn planar_gaussian_f_is_closed_form_and_fg_is_additive() {
    let (mu, sigma) = (2.0, 1.3);
    // Unit normal along +x, plane through origin ⇒ ξ = x-coordinate.
    let coord = Coordinate::planar([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let target = ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian { mu, sigma }));
    let r = ProfileRestraint::new(coord, target, KT, PbcWrap::none());

    let scale = 0.75;
    let k_u = KT / (2.0 * sigma * sigma); // U = k_u·(ξ−μ)²
    let k_du = KT / (sigma * sigma); // dU/dξ = k_du·(ξ−μ)

    for &xq in &[0.5, 1.0, 2.0, 3.5, 5.0] {
        let x = [xq, 7.0, -3.0]; // y, z are irrelevant to a +x planar coordinate
        let xi = xq; // ξ = x-coordinate

        // f closed form: scale·(kT/2σ²)(ξ−μ)².
        let want_f = scale * k_u * (xi - mu).powi(2);
        let got_f = r.f(&x, scale, 0.0);
        assert!(
            (got_f - want_f).abs() <= TOL * want_f.abs().max(1.0),
            "f @ {xq}: got {got_f} want {want_f}"
        );

        // fg returns the same value AND accumulates into a PRE-SEEDED g.
        let seed = [0.11, -0.22, 0.33];
        let mut g = seed;
        let got_fg = r.fg(&x, scale, 0.0, &mut g);
        assert!(
            (got_fg - want_f).abs() <= TOL * want_f.abs().max(1.0),
            "fg value @ {xq}: got {got_fg} want {want_f}"
        );

        // n̂ = (1,0,0) ⇒ gradient adds only on x; additive (seed + contribution).
        let contrib_x = scale * k_du * (xi - mu);
        assert!(
            (g[0] - (seed[0] + contrib_x)).abs() <= TOL,
            "g[0] not additive @ {xq}: got {} want {}",
            g[0],
            seed[0] + contrib_x
        );
        assert!((g[1] - seed[1]).abs() <= TOL, "g[1] overwritten @ {xq}");
        assert!((g[2] - seed[2]).abs() <= TOL, "g[2] overwritten @ {xq}");
    }
}

// ── ac-002: full 4×5 coordinate×distribution gradient matrix vs FD ──────────

/// Central finite difference of `f` componentwise.
fn fd_grad(r: &ProfileRestraint, x: &[F; 3], scale: F, h: F) -> [F; 3] {
    let mut g = [0.0; 3];
    for k in 0..3 {
        let mut xp = *x;
        xp[k] += h;
        let mut xm = *x;
        xm[k] -= h;
        g[k] = (r.f(&xp, scale, 0.0) - r.f(&xm, scale, 0.0)) / (2.0 * h);
    }
    g
}

/// Assert the analytic `fg` gradient agrees with central FD of `f` componentwise
/// at a sample point chosen away from singularities and the density floor.
fn assert_fg_matches_fd(r: &ProfileRestraint, x: &[F; 3], scale: F, label: &str) {
    let analytic = {
        let mut g = [0.0; 3];
        r.fg(x, scale, 0.0, &mut g);
        g
    };
    let fd = fd_grad(r, x, scale, 1e-6);
    for k in 0..3 {
        let denom = analytic[k].abs().max(fd[k].abs()).max(1.0);
        let rel = (analytic[k] - fd[k]).abs() / denom;
        assert!(
            rel <= 1e-4,
            "{label}: axis {k} analytic={} fd={} rel_err={}",
            analytic[k],
            fd[k],
            rel
        );
    }
}

/// The five distribution targets for the matrix, parameterised so the chosen
/// sample point lands in a smooth, well-conditioned region for each coordinate.
fn distribution_targets(mu: F) -> Vec<(&'static str, ProfileTarget)> {
    vec![
        (
            "gaussian",
            ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian { mu, sigma: 1.5 })),
        ),
        (
            "erf",
            ProfileTarget::Analytic(planar_penalty(Distribution::ErfStep {
                xi0: mu,
                w: 1.2,
                rising: true,
            })),
        ),
        (
            "tanh",
            ProfileTarget::Analytic(planar_penalty(Distribution::TanhStep {
                xi0: mu,
                w: 1.2,
                rising: true,
            })),
        ),
        (
            "exponential",
            ProfileTarget::Analytic(planar_penalty(Distribution::Exponential { lambda: 2.5 })),
        ),
        (
            "tabulated",
            ProfileTarget::Tabulated(tabulated_gaussian(mu, 1.5)),
        ),
    ]
}

#[test]
fn gradient_matrix_matches_finite_difference() {
    let scale = 0.9;
    let pbc = PbcWrap::none();

    // Each coordinate paired with a sample point whose ξ is positive, off the
    // singularity, and inside the tabulated grid (mu±4σ around mu=5 ⇒ ξ∈[-1,11]).
    let mu = 5.0;

    // planar: ξ = x-coordinate (normal +x, plane through origin).
    let planar = Coordinate::planar([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let planar_pt = [5.5, 2.0, -1.0]; // ξ = 5.5

    // radial: centre at origin; point at radius 5 (well outside r_guard).
    let radial = Coordinate::radial([0.0; 3], 0.1);
    let radial_pt = [3.0, 4.0, 0.0]; // ξ = 5.0

    // cylindrical: z-axis through origin; perpendicular distance 5.
    let cylindrical = Coordinate::cylindrical([0.0; 3], [0.0, 0.0, 1.0], 0.1).unwrap();
    let cyl_pt = [3.0, 4.0, 8.0]; // ξ = sqrt(9+16) = 5

    // region-distance: box [0,10]³; point outside +x face ⇒ signed distance 5.
    let region = Arc::new(InsideBoxRegion::new([0.0; 3], [10.0; 3]));
    let region_coord = Coordinate::region_distance(region);
    let region_pt = [15.0, 5.0, 5.0]; // signed distance = 5 (smooth, single face)

    let coords: Vec<(&str, Coordinate, [F; 3])> = vec![
        ("planar", planar, planar_pt),
        ("radial", radial, radial_pt),
        ("cylindrical", cylindrical, cyl_pt),
        ("region", region_coord, region_pt),
    ];

    for (cname, coord, pt) in &coords {
        for (dname, target) in distribution_targets(mu) {
            let r = ProfileRestraint::new(coord.clone(), target, KT, pbc);
            assert_fg_matches_fd(&r, pt, scale, &format!("{cname}×{dname}"));
        }
    }
}

// ── ac-003: guard + floor keep f / g finite at singular sites ───────────────

#[test]
fn finite_at_radial_centre_and_cylindrical_axis() {
    let target = ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian {
        mu: 3.0,
        sigma: 1.5,
    }));

    // Radial centre exactly on the singularity (inside r_guard).
    let radial = Coordinate::radial([1.0, 2.0, 3.0], 0.25);
    let r_rad = ProfileRestraint::new(radial, target.clone(), KT, PbcWrap::none());
    let centre = [1.0, 2.0, 3.0];
    assert_all_finite(&r_rad, &centre, "radial centre");

    // Cylindrical axis exactly on the singularity (inside r_guard).
    let cyl = Coordinate::cylindrical([0.0; 3], [0.0, 0.0, 1.0], 0.25).unwrap();
    let r_cyl = ProfileRestraint::new(cyl, target, KT, PbcWrap::none());
    let on_axis = [0.0, 0.0, 7.0];
    assert_all_finite(&r_cyl, &on_axis, "cylindrical axis");
}

#[test]
fn finite_in_zero_density_bin() {
    // Tabulated profile with empty middle bins ⇒ ξ in those bins is floored.
    let nodes = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let values = vec![1.0, 0.0, 0.0, 0.0, 1.0]; // empty interior bins
    let tab = TabulatedProfile::new(
        nodes,
        values,
        InputKind::VolumetricDensity,
        ShellJacobian::Planar,
        DensityFloor::new(1e-3, 1.0, KT),
    )
    .expect("valid grid");
    let coord = Coordinate::planar([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let r = ProfileRestraint::new(coord, ProfileTarget::Tabulated(tab), KT, PbcWrap::none());

    // ξ = 2.0 lands in the floored zero-density region.
    let x = [2.0, 0.0, 0.0];
    assert_all_finite(&r, &x, "zero-density bin");
}

/// Assert `f` and every component of `g` (from `fg`) are finite (no Inf/NaN).
fn assert_all_finite(r: &ProfileRestraint, x: &[F; 3], label: &str) {
    let f = r.f(x, 1.0, 0.0);
    assert!(f.is_finite(), "{label}: f not finite: {f}");
    let mut g = [0.0; 3];
    let fg = r.fg(x, 1.0, 0.0, &mut g);
    assert!(fg.is_finite(), "{label}: fg value not finite: {fg}");
    for (k, gk) in g.iter().enumerate() {
        assert!(gk.is_finite(), "{label}: g[{k}] not finite: {gk}");
    }
}

// ── ac-004: Send + Sync + Debug; linear in scale; periodic_box mirrors box ──

#[test]
fn is_send_sync_and_debug() {
    fn _assert_send_sync<T: Send + Sync>() {}
    _assert_send_sync::<ProfileRestraint>();
    let coord = Coordinate::planar([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let target = ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian {
        mu: 1.0,
        sigma: 1.0,
    }));
    let r = ProfileRestraint::new(coord, target, KT, PbcWrap::none());
    let _ = format!("{r:?}");
}

#[test]
fn doubling_scale_doubles_f_and_g() {
    let coord = Coordinate::radial([0.0; 3], 0.1);
    let target = ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian {
        mu: 3.0,
        sigma: 1.5,
    }));
    let r = ProfileRestraint::new(coord, target, KT, PbcWrap::none());
    let x = [3.0, 4.0, 0.0]; // ξ = 5

    let f1 = r.f(&x, 1.0, 0.0);
    let f2 = r.f(&x, 2.0, 0.0);
    assert!(
        (f2 - 2.0 * f1).abs() <= TOL * f1.abs().max(1.0),
        "f not linear in scale: f1={f1} f2={f2}"
    );

    // Accumulated g must double too (with a zero seed, the contribution doubles).
    let mut g1 = [0.0; 3];
    r.fg(&x, 1.0, 0.0, &mut g1);
    let mut g2 = [0.0; 3];
    r.fg(&x, 2.0, 0.0, &mut g2);
    for k in 0..3 {
        assert!(
            (g2[k] - 2.0 * g1[k]).abs() <= TOL * g1[k].abs().max(1.0),
            "g[{k}] not linear in scale: g1={} g2={}",
            g1[k],
            g2[k]
        );
    }
}

#[test]
fn periodic_box_mirrors_pbc_when_anchored() {
    let coord = Coordinate::planar([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
    let target = ProfileTarget::Analytic(planar_penalty(Distribution::Gaussian {
        mu: 1.0,
        sigma: 1.0,
    }));

    // Non-periodic carrier ⇒ no box.
    let r_none = ProfileRestraint::new(coord.clone(), target.clone(), KT, PbcWrap::none());
    assert!(r_none.periodic_box().is_none(), "non-periodic must be None");

    // Periodic on x only ⇒ Some box mirroring the pbc lengths and flags.
    let pbc = PbcWrap::new([10.0, 12.0, 14.0], [true, false, false]);
    let r_pbc = ProfileRestraint::new(coord, target, KT, pbc);
    let (min, max, periodic) = r_pbc.periodic_box().expect("periodic axis ⇒ Some");
    assert_eq!(min, [0.0, 0.0, 0.0]);
    assert_eq!(max, [10.0, 12.0, 14.0]);
    assert_eq!(periodic, [true, false, false]);
}
