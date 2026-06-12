---
slug: profile-distribution-restraints-03-distribution
criteria:
  - id: ac-001
    summary: Closed-form du_dxi matches central FD of u for all four shapes
    type: code
    evaluator_hint: "test module in src/restraint/profile/distribution.rs"
    pass_when: |
      A test in src/restraint/profile/distribution.rs central-differences
      u(xi, kt) at representative xi and asserts |du_dxi − fd| within the
      stated tolerance (e.g. rel <= 1e-6 at h≈1e-5) for Gaussian, ErfStep,
      TanhStep, and Exponential; the test passes.
    status: pending
  - id: ac-002
    summary: Gaussian equals harmonic forms; exponential force is constant kT/λ
    type: code
    evaluator_hint: "test module in src/restraint/profile/distribution.rs"
    pass_when: |
      A test asserts Gaussian u == (kT/2σ²)(ξ−μ)² and du_dxi == (kT/σ²)(ξ−μ)
      to floating tolerance, and Exponential du_dxi == kT/λ (constant,
      independent of ξ) at several ξ≥0; the test passes.
    status: pending
  - id: ac-003
    summary: Radial histogram-with-Jacobian equals density U; omitting it differs
    type: code
    evaluator_hint: "test module in src/restraint/profile/distribution.rs"
    pass_when: |
      A two-sided test builds n(ξ) and ρ*(ξ)=n(ξ)/(4π ξ²) for the same radial
      profile: ProfilePenalty(CountHistogram, Radial) U(ξ) equals
      ProfilePenalty(VolumetricDensity) U(ξ) to floating tolerance, AND
      treating the histogram as a density yields a different U whose dU/dξ
      differs by the predicted +2kT/ξ; the test passes.
    status: pending
  - id: ac-004
    summary: Zero-density region gives finite capped U and finite du_dxi
    type: code
    evaluator_hint: "test module in src/restraint/profile/distribution.rs"
    pass_when: |
      A test drives the target density below ρ_min (or to a numerically-zero
      erf/tanh tail) and asserts U is finite and equals
      U_max = −kT·ln(ρ_min/ρ₀), and du_dxi is finite with no Inf/NaN; the
      test passes.
    status: pending
---

# Acceptance criteria

- **ac-001 — analytic gradient is correct.** The whole point of shipping closed forms instead of numerically differentiating is that the analytic `du_dxi` must agree with a finite-difference of `u`. Verified per shape at representative points away from the floor-capped region.
- **ac-002 — known limits hold exactly.** Gaussian must reduce to the harmonic well (spring `k = kT/σ²` centred at μ) and the exponential must give a constant force `kT/λ`; these pin the constants and catch sign/factor errors a noisy FD check could miss.
- **ac-003 — the load-bearing §6.2 Jacobian.** The histogram→density conversion is the single most failure-prone correctness point: a count histogram and the matching volumetric density must invert to the same `U` only when the shell factor (∝ξ² radial) is divided out, and omitting it must produce the predicted spurious `+2kT/ξ` outward force. Tested two-sided so a no-op implementation cannot pass.
- **ac-004 — divergence control.** A zero-probability bin must be capped at `U_max = −kT·ln(ρ_min/ρ₀)` with a finite gradient, never an `Inf`/`NaN` wall that would stall gencan.
