//! Torsion priors, calibrated against the chain statistics they exist to
//! reproduce.

use super::*;

/// `three_state_from_c_inf(5.5, 109.47°)` — the PEO calibration (spec §5.1):
/// States with trans at |φ| = π and gauche± at ±π/3 (IUPAC absolute
/// convention), weights (p_t, p_g, p_g) normalized, p_t ≈ 0.645 solved from
/// C∞ = C_FRC·(1+⟨cosφ'⟩)/(1−⟨cosφ'⟩) with φ' measured from trans.
#[test]
fn prior_three_state_calibration() {
    let pi = std::f64::consts::PI as F;
    let theta = (109.47 as F).to_radians();
    let prior = TorsionPrior::three_state_from_c_inf(5.5, theta);
    let TorsionPrior::States(states) = &prior else {
        panic!("three_state_from_c_inf must return TorsionPrior::States, got {prior:?}");
    };
    assert_eq!(states.len(), 3, "three states: trans + gauche±");
    let total: F = states.iter().map(|&(_, w)| w).sum();
    assert!(
        (total - 1.0).abs() < 1e-12,
        "weights must sum to 1, got {total}"
    );

    let (trans, mut gauche): (Vec<_>, Vec<_>) = states
        .iter()
        .copied()
        .partition(|&(a, _)| (a.abs() - pi).abs() < 1e-9);
    assert_eq!(
        trans.len(),
        1,
        "exactly one trans state at |φ| = π, got {states:?}"
    );
    // C_FRC = (1−cosθ)/(1+cosθ) ≈ 2.000; r = 5.5/C_FRC ≈ 2.750;
    // x = (r−1)/(r+1) ≈ 0.4667; p_t = (2x+1)/3 ≈ 0.6445.
    assert!(
        (trans[0].1 - 0.645).abs() < 0.01,
        "trans weight = {}, expected ≈ 0.645 for C∞ = 5.5",
        trans[0].1
    );
    gauche.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("finite angles"));
    assert_eq!(gauche.len(), 2, "two gauche states");
    assert!(
        (gauche[0].0 + pi / 3.0).abs() < 1e-9,
        "gauche− at −π/3, got {}",
        gauche[0].0
    );
    assert!(
        (gauche[1].0 - pi / 3.0).abs() < 1e-9,
        "gauche+ at +π/3, got {}",
        gauche[1].0
    );
    assert!(
        (gauche[0].1 - gauche[1].1).abs() < 1e-12,
        "gauche weights must be equal: {} vs {}",
        gauche[0].1,
        gauche[1].1
    );
}

/// The ac-006(1) regression baseline: uniform torsions on a fixed-109.5°
/// chain are the freely rotating chain, C∞ = (1−cosθ)/(1+cosθ) = 2.00
/// exactly. This pins the NeRF rebuild + sampler combination to an analytic
/// value; the fixed seed makes the sampled number stable.
#[test]
fn prior_uniform_freely_rotating_c_inf() {
    let c_n = sampled_c_n(&TorsionPrior::Uniform, 200, 1.53, 600, 42);
    assert!(
        (c_n - 2.0).abs() <= 0.1,
        "freely rotating chain: C_n = {c_n:.4}, analytic C∞ = 2.00 ± 0.1 (spec §5.1)"
    );
}

/// `States([(π, 1)])` is the all-trans prior, and the zigzag template IS the
/// all-trans conformer — sampling must reproduce the template coordinates.
#[test]
fn prior_states_trans_only_is_all_trans() {
    let pi = std::f64::consts::PI as F;
    let template = zigzag_coords(20, 1.53);
    let tree = InternalTree::from_frame(
        &frame_from_parts(&template, &chain_bonds(20)),
        &BondDistanceWeights::from_exclusion_depth(3),
    )
    .expect("chain template decomposes");
    let prior = TorsionPrior::States(vec![(pi, 1.0)]);
    let mut rng = SmallRng::seed_from_u64(7);
    let mut vars = vec![0.0 as F; tree.n_vars()];
    for k in 0..tree.n_steps() {
        if let Some(v) = tree.step_var(k) {
            vars[v] = prior.sample(tree.template_var(k), &mut rng);
        }
    }
    let rebuilt = rebuild_coords(&tree, &template, &vars);
    let dev = max_abs_dev(&template, &rebuilt);
    assert!(
        dev < 1e-6,
        "an all-trans States prior must reproduce the (already all-trans) \
         zigzag template: ‖Δ‖∞ = {dev:e}"
    );
}

/// The calibration round-trip (spec chain-statistics ladder, step 3): a chain
/// sampled from `three_state_from_c_inf(5.5, 109.47°)` must measure back
/// C_n = 5.5 ± 0.3. Finite-n droop is real at n = 200; if this fails
/// marginally when Task 4 lands, raise n_beads to 400 and keep the tolerance.
#[test]
fn prior_ris_calibrated_c_inf() {
    let theta = (109.47 as F).to_radians();
    let prior = TorsionPrior::three_state_from_c_inf(5.5, theta);
    let c_n = sampled_c_n(&prior, 200, 1.53, 600, 42);
    assert!(
        (c_n - 5.5).abs() <= 0.3,
        "RIS-calibrated chain: C_n = {c_n:.4}, target C∞ = 5.5 ± 0.3"
    );
}
