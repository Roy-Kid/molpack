---
slug: profile-distribution-restraints-06-script
criteria:
  - id: ac-001
    summary: profile keyword parses to RestraintSpec::Profile with atoms-subset scoping
    type: code
    evaluator_hint: "test: parser round-trip in src/script/parser.rs"
    pass_when: |
      A round-trip unit test parses a `.inp` `profile gaussian plane <n> <x0>
      mu <μ> sigma <σ>` line into `RestraintSpec::Profile` with the matching
      distribution, geometry, and params and `input_kind == Histogram` by
      default; the same line placed inside `atoms <idx…> … end atoms` lands on
      that AtomGroup's restraints scoped to the listed indices, and the test
      passes under `cargo test -p molcrafts-molpack --lib`.
    status: pending
  - id: ac-002
    summary: lowered profile spec yields a ProfileRestraint matching keyword args
    type: code
    evaluator_hint: "test: lowering check via restraint_from_spec"
    pass_when: |
      A test lowers a parsed `RestraintSpec::Profile` through
      `restraint_from_spec` and asserts the resulting `ProfileRestraint`'s
      coordinate and distribution reproduce the keyword arguments — verified by
      its `f`/`fg` output at known ξ points matching the analytic form within
      1e-9 — and passes under `cargo test -p molcrafts-molpack --lib --tests`.
    status: pending
  - id: ac-003
    summary: planar Gaussian packs in tolerance; radial histogram needs the ξ² Jacobian
    type: code
    evaluator_hint: "test: profile_pack_e2e.rs measured-profile helper"
    pass_when: |
      In `tests/profile_pack_e2e.rs`, a single-component planar-Gaussian `.inp`
      packs to a measured ρ(z) within the test's stated tolerance, and a radial
      count-histogram `.inp` matches its target within tolerance ONLY with the
      shell-volume (∝ξ²) Jacobian applied and is asserted out of tolerance
      without it (two-sided); both assertions pass under
      `cargo test -p molcrafts-molpack --tests`.
    status: pending
  - id: ac-004
    summary: opposite-erf gives asymmetric layout; zero-density region packs without NaN
    type: code
    evaluator_hint: "test: profile_pack_e2e.rs two-component + zero-density cases"
    pass_when: |
      In `tests/profile_pack_e2e.rs`, a two-component opposite-erf `.inp` yields
      per-component measured profiles each biased to its own side (asymmetric
      layout asserted), and a target containing a zero-density region packs to
      completion with all final positions and gradients finite (no NaN);
      both pass under `cargo test -p molcrafts-molpack --tests`.
    status: pending
---

# Acceptance criteria

- **ac-001 — parser round-trip.** Proves the new `profile` keyword is wired into the script grammar at both scopes (`State::InStructure` and `State::InAtoms`) and captures the density-vs-histogram flag with its default. Covers Tasks 1, 4, 5.
- **ac-002 — lowering fidelity.** Proves the `restraint_from_spec` arm reconstructs the correct `-05` coordinate+distribution composition, validated against the analytic `f`/`fg`. Covers Tasks 2, 6.
- **ac-003 — planar + radial domain validation.** The first two user-facing packing results; the radial case is the two-sided §6.2 check that the Jacobian is both present and load-bearing. Covers Task 3.
- **ac-004 — interface superposition + numerical robustness.** The §6.7 two-component asymmetry plus the zero-density NaN-safety guarantee. Covers Task 3.
