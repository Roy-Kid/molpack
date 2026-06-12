---
slug: profile-distribution-restraints-04-spline
criteria:
  - id: ac-001
    summary: Spline interpolates nodes exactly and is C1 across interior nodes
    type: code
    evaluator_hint: "cargo test -p molcrafts-molpack --lib spline::"
    pass_when: |
      Tests in src/restraint/profile/spline.rs assert TabulatedProfile::value
      returns each tabulated rho*_i at its node xi_i within f64 round-off, and
      that value and deriv from the left- and right-adjacent segments agree at
      every interior node within tolerance (C1 continuity). Tests pass.
    status: verified
  - id: ac-002
    summary: Analytic spline derivative matches central finite difference
    type: code
    evaluator_hint: "cargo test -p molcrafts-molpack --lib spline::"
    pass_when: |
      A test compares TabulatedProfile::deriv(xi) against a central finite
      difference of TabulatedProfile::value across the grid interior and asserts
      agreement within the stated tolerance. Test passes.
    status: verified
  - id: ac-003
    summary: Tabulated Gaussian reproduces analytic harmonic U within tolerance
    type: code
    evaluator_hint: "cargo test -p molcrafts-molpack --lib spline::"
    pass_when: |
      A test samples a Gaussian density rho* ∝ exp(-(xi-mu)^2/2sigma^2) onto a
      grid, builds TabulatedProfile, and asserts u(xi, kt) matches the analytic
      harmonic U(xi) = (kT/2sigma^2)(xi-mu)^2 (up to the additive rho0 constant)
      within the stated interpolation tolerance. Test passes.
    status: verified
  - id: ac-004
    summary: Monotone input stays overshoot-free; floor and radial Jacobian hold
    type: code
    evaluator_hint: "cargo test -p molcrafts-molpack --lib spline::"
    pass_when: |
      Tests assert (a) a monotone-decreasing tabulated input yields a spline
      whose values never leave the local node bracket (no overshoot), (b) a
      zero-valued bin yields finite u and finite du_dxi via the density floor,
      and (c) a radial count histogram reproduces the intended uniform density
      only when the 4*pi*xi^2 shell Jacobian is applied. All tests pass.
    status: verified
---

# Acceptance criteria

- **ac-001** — Exact node interpolation plus C¹ continuity. The spline must pass through every `(ξ_i, ρ*_i)`, and both value and first derivative must be continuous across each interior node. This is the structural correctness bar for the interpolant.
- **ac-002** — The stored per-segment analytic derivative `s'(ξ)` must match a central finite difference of `s(ξ)` across the grid interior within a stated tolerance, proving `deriv` is the true derivative of `value` (required so `du_dxi` is a faithful force).
- **ac-003** — Special-case identity: a tabulated sampling of a Gaussian density, after Boltzmann inversion, reproduces the analytic harmonic `U(ξ) = (kT/2σ²)(ξ−μ)²` (modulo the ρ₀ additive constant) within interpolation tolerance. Confirms the analytic families are special cases of the general tabulated profile.
- **ac-004** — Robustness of the inversion: monotone input produces an overshoot-free monotone spline (Fritsch–Carlson limiter working); a zero bin yields finite `u` and `du_dxi` (density floor working); and a radial histogram only reproduces the intended density when the ∝ξ² shell Jacobian is applied (§6.2 conversion working).
