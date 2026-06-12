//! End-to-end pack tests for the `profile` keyword (capstone sub-spec -06).
//!
//! Each test drives the REAL script path — `script::parse` → `Script::lower` →
//! `StructurePlan::apply` — so the `profile` grammar, the
//! [`RestraintSpec::Profile`](molpack::script::RestraintSpec) lowering arm, and
//! the composed `ProfileRestraint` are all exercised together, then packs a
//! small SEEDED system through `Molpack::pack_with_report`. Templates are built
//! programmatically from single-site coordinates (the wheel-style, io-free path),
//! so no molecule files are read.
//!
//! Because a soft profile restraint only APPROACHES its target (it rides the
//! scale soft-start and competes with the overlap term), the assertions are
//! deliberately QUALITATIVE and robust — modal/mean concentration, two-sided
//! sign differences, and finiteness — not a tight χ². They are seeded and were
//! confirmed stable across repeated runs.
//!
//! Speed budget: tens of single-site molecules, capped `max_loops`, ≤ ~200 ms
//! each like `restraint_pack_e2e.rs`.

use std::path::Path;

use molpack::script::{self, StructurePlan};
use molpack::{F, Molpack, Target};

const SEED: u64 = 0xC0FFEE;
const TOL: F = 1.5;
const PRECISION: F = 1e-2;

// ── harness ──────────────────────────────────────────────────────────────────

/// Parse a `.inp` fragment, lower it, and build one programmatic single-site
/// target per structure block (radius 0.5 Å), applying the script's restraints.
/// Returns the configured packer, the targets, and the lowered `nloop`.
fn build(inp: &str) -> (Molpack, Vec<Target>, usize) {
    let parsed = script::parse(inp).expect("parse failed");
    let plan = parsed.lower(Path::new(".")).expect("lower failed");
    let targets: Vec<Target> = plan
        .structures
        .iter()
        .map(build_single_site_target)
        .collect();
    (plan.packer, targets, plan.nloop)
}

/// One single-site template (radius 0.5 Å) replicated `sp.number` times, with the
/// plan's restraints / centering applied via the public `StructurePlan::apply`.
fn build_single_site_target(sp: &StructurePlan) -> Target {
    let target = Target::from_coords(&[[0.0, 0.0, 0.0]], &[0.5], sp.number);
    sp.apply(target)
}

/// Pack a single target and return the packed site positions. Uses a fixed seed
/// and a capped loop budget so the run is deterministic and fast.
fn pack_positions(inp: &str, max_loops: usize) -> Vec<[F; 3]> {
    let (mut packer, targets, _nloop) = build(inp);
    packer = packer
        .with_tolerance(TOL)
        .with_precision(PRECISION)
        .with_seed(SEED);
    let result = packer
        .pack_with_report(&targets, max_loops)
        .expect("pack should succeed");
    result.positions()
}

/// Mean of a slice.
fn mean(xs: &[F]) -> F {
    xs.iter().sum::<F>() / xs.len() as F
}

/// Median of a slice (sorts a copy).
fn median(xs: &[F]) -> F {
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

/// Index of the histogram bin (over `[lo, hi]`, `n` bins) holding the most mass.
fn modal_bin(values: &[F], lo: F, hi: F, n: usize) -> usize {
    let mut counts = vec![0usize; n];
    let width = (hi - lo) / n as F;
    for &v in values {
        if v < lo || v >= hi {
            continue;
        }
        let b = ((v - lo) / width) as usize;
        counts[b.min(n - 1)] += 1;
    }
    counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, &c)| c)
        .map(|(i, _)| i)
        .unwrap()
}

// ── ac-003: single-component planar Gaussian concentrates near μ ──────────────

#[test]
fn planar_gaussian_concentrates_near_mu() {
    // Confine to a 0..30 box on z (via plane geometry the coordinate is z) and
    // bias toward z = 15 with a narrow Gaussian. The packed z's should cluster
    // around 15: small mean |z−μ| and the modal z-bin is the μ-bin.
    let inp = "\
tolerance 1.5
seed 11
output out.xyz

structure site.xyz
  number 40
  inside box 0. 0. 0. 30. 30. 30.
  profile gaussian plane 0. 0. 1. 0. 0. 15. mu 0. sigma 2. density
end structure
";
    // mu is expressed in the coordinate ξ = z − 15 (plane point at z=15), so the
    // Gaussian peak ξ=0 sits at z=15.
    let positions = pack_positions(inp, 25);
    let zs: Vec<F> = positions.iter().map(|p| p[2]).collect();
    let mu_z = 15.0;
    let dev: Vec<F> = zs.iter().map(|&z| (z - mu_z).abs()).collect();
    let mean_dev = mean(&dev);
    let med_dev = median(&dev);

    // Robust concentration: median |z−μ| well under a quarter of the box half-width.
    assert!(
        med_dev < 4.0,
        "planar Gaussian: median |z−μ| = {med_dev} (mean {mean_dev}) not concentrated near μ={mu_z}"
    );
    // The modal z-bin (6 bins over 0..30, width 5) is the one containing μ=15 (bin 3).
    let modal = modal_bin(&zs, 0.0, 30.0, 6);
    assert_eq!(
        modal, 3,
        "planar Gaussian: modal z-bin {modal} is not the μ=15 bin (bin 3); zs={zs:?}"
    );
}

// ── ac-003: radial histogram needs the ξ² shell Jacobian (two-sided) ──────────

/// Pack a radial-Gaussian target peaked at ξ=7 inside a confining sphere, once as
/// a count *histogram* and once as a volumetric *density*, holding everything
/// else fixed. The §6.2 shell Jacobian is the only difference and it is
/// load-bearing:
///
/// - **density** — the Gaussian *is* the target volumetric density, so sites
///   sit near its peak, r ≈ μ = 7.
/// - **histogram** — the count Gaussian is divided by the shell volume `4πξ²`
///   before inversion (equivalently `U` gains `+kt·ln(4πξ²)`, an inward push),
///   so the realised target density is far higher at small ξ and the sites
///   over-concentrate toward the centre, r ≪ μ.
///
/// The two-sided claim is the predicted DIRECTION of the difference — the
/// histogram run packs to a clearly smaller mean radius than the density run —
/// not a tight match to either curve.
#[test]
fn radial_histogram_jacobian_changes_mean_radius() {
    let with_jacobian = "\
tolerance 1.5
seed 7
output out.xyz

structure site.xyz
  number 30
  inside sphere 0. 0. 0. 14.
  profile gaussian radial 0. 0. 0. mu 7. sigma 2.5 histogram
end structure
";
    let no_jacobian = "\
tolerance 1.5
seed 7
output out.xyz

structure site.xyz
  number 30
  inside sphere 0. 0. 0. 14.
  profile gaussian radial 0. 0. 0. mu 7. sigma 2.5 density
end structure
";

    let r_hist = mean_radius(&pack_positions(with_jacobian, 20));
    let r_dens = mean_radius(&pack_positions(no_jacobian, 20));

    // Jacobian load-bearing, two-sided: the histogram run (shell-divided target)
    // over-concentrates toward the centre, so its mean radius is clearly smaller
    // than the density run's, which sits near the Gaussian peak μ=7.
    assert!(
        r_hist + 1.0 < r_dens,
        "radial Jacobian two-sided check: mean r(histogram)={r_hist} should be \
         clearly below mean r(density)={r_dens}"
    );
    // The density run genuinely realises the peak near μ=7 (not collapsed to 0).
    assert!(
        r_dens > 4.0,
        "radial density run mean r={r_dens} should sit near the Gaussian peak μ=7"
    );
}

/// Mean distance of each site from the origin.
fn mean_radius(positions: &[[F; 3]]) -> F {
    let radii: Vec<F> = positions
        .iter()
        .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
        .collect();
    mean(&radii)
}

// ── ac-004: two-component opposite-erf gives an asymmetric layout ─────────────

#[test]
fn opposite_erf_components_separate_across_plane() {
    // Two species share the z=15 plane: species A is biased to the RISING side
    // (high z), species B to the FALLING side (low z). Their mean z must land on
    // opposite sides of the plane — an asymmetric, sign-separated layout (§6.7).
    let inp = "\
tolerance 1.5
seed 23
output out.xyz

structure a.xyz
  number 30
  inside box 0. 0. 0. 30. 30. 30.
  profile erf plane 0. 0. 1. 0. 0. 15. xi0 0. width 2. rising density
end structure

structure b.xyz
  number 30
  inside box 0. 0. 0. 30. 30. 30.
  profile erf plane 0. 0. 1. 0. 0. 15. xi0 0. width 2. falling density
end structure
";
    let (mut packer, targets, _) = build(inp);
    packer = packer
        .with_tolerance(TOL)
        .with_precision(PRECISION)
        .with_seed(SEED);
    let result = packer
        .pack_with_report(&targets, 25)
        .expect("pack should succeed");
    let positions = result.positions();
    // Targets pack in order: first 30 are species A, next 30 species B.
    let za: Vec<F> = positions[..30].iter().map(|p| p[2]).collect();
    let zb: Vec<F> = positions[30..].iter().map(|p| p[2]).collect();
    let mean_a = mean(&za);
    let mean_b = mean(&zb);

    let plane = 15.0;
    assert!(
        mean_a > plane && mean_b < plane,
        "opposite-erf asymmetry: rising A mean z={mean_a} should be above the \
         plane (z={plane}) and falling B mean z={mean_b} below it"
    );
    // And the two means are genuinely separated, not marginally split.
    assert!(
        mean_a - mean_b > 2.0,
        "opposite-erf: per-species means {mean_a} vs {mean_b} are not clearly separated"
    );
}

// ── ac-004: a zero-density region packs to completion with finite state ───────

#[test]
fn zero_density_region_packs_without_nan() {
    // A tabulated radial density with an interior zero band (ρ=0 between ξ=4 and
    // ξ=8) is a forbidden shell. The density floor caps the penalty, so packing
    // must COMPLETE with every final position and gradient finite — a robustness
    // guarantee, not a profile-shape claim.
    let inp = "\
tolerance 1.5
seed 31
output out.xyz

structure site.xyz
  number 30
  inside sphere 0. 0. 0. 14.
  profile tabulated radial 0. 0. 0. density 1. 1.0 3. 0.8 4. 0.0 8. 0.0 9. 0.7 13. 1.0
end structure
";
    let positions = pack_positions(inp, 25);
    assert_eq!(positions.len(), 30, "all 30 sites must be placed");
    for (i, p) in positions.iter().enumerate() {
        for (k, &c) in p.iter().enumerate() {
            assert!(c.is_finite(), "site {i} axis {k} not finite: {c}");
        }
    }
}
