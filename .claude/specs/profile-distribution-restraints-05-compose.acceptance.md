---
slug: profile-distribution-restraints-05-compose
criteria:
  - id: ac-001
    summary: planar-Gaussian f matches closed form; fg accumulates additively
    type: code
    evaluator_hint: "module: restraint::profile::tests::compose"
    pass_when: |
      A test in src/restraint/profile/tests/compose.rs asserts ProfileRestraint::f
      with a planar Gaussian returns scale·(kT/2σ²)(ξ−μ)² at sample points, and fg
      returns that same value while leaving g equal to its pre-seeded contents plus
      scale·(kT/σ²)(ξ−μ)·n̂ (additive, never overwritten); test passes.
    status: verified
  - id: ac-002
    summary: full coordinate×distribution analytic gradient matches finite difference
    type: code
    evaluator_hint: "module: restraint::profile::tests::compose"
    pass_when: |
      A parameterised test covers every pair in {planar, radial, cylindrical,
      region-distance} × {Gaussian, erf, tanh, exponential, tabulated} and asserts
      the fg analytic gradient agrees with a central finite difference of f
      componentwise within the stated tolerance for all 20 pairs; test passes.
    status: verified
  - id: ac-003
    summary: guard + floor compose to keep f and g finite at singular sites
    type: code
    evaluator_hint: "module: restraint::profile::tests::compose"
    pass_when: |
      A test places a biased site exactly on the radial centre / cylindrical axis
      (inside r_guard) and another in a zero-density bin, and asserts f and every
      component of g are finite (no Inf/NaN) in both cases; test passes.
    status: verified
  - id: ac-004
    summary: Send+Sync+Debug, no interior mutability, linear in global scale
    type: code
    evaluator_hint: "module: restraint::profile::tests::compose"
    pass_when: |
      ProfileRestraint compiles under a static assert_send_sync bound and derives
      Debug; a test asserts that doubling scale exactly doubles both f and the
      accumulated g (and periodic_box mirrors the coordinate's anchored box);
      test passes and the type holds no Cell/RefCell/Mutex/atomic field.
    status: verified
---

# Acceptance criteria

- **ac-001** — Pins the composition's value/gradient correctness on the simplest analytic case and the load-bearing additivity contract (`g += …`, never overwrite), matching the objective's accumulation loop at src/objective.rs:544–546.
- **ac-002** — The full §6.9 gradient matrix; analytic-vs-FD agreement across all 4×5 coordinate×distribution pairs is the binding proof that the chain rule ∇ₓU = (dU/dξ)·∇ξ is composed correctly for every member of the family.
- **ac-003** — Confirms the -02 inner clamp and the -03/-04 floor compose so a site on a coordinate singularity or in an empty density bin never produces a non-finite force.
- **ac-004** — Locks the soft-start coupling (§6.8): linear scaling through `scale` with no mutable ramp state, plus the `Send + Sync + Debug` requirements for parallel evaluation and periodic-box consistency.
