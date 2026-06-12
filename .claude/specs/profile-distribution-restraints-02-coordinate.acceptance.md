---
slug: profile-distribution-restraints-02-coordinate
criteria:
  - id: ac-001
    summary: Analytic grad_xi matches central FD of xi componentwise for all four coordinates
    type: code
    evaluator_hint: "test target: tests/coordinate.rs analytic-vs-FD cases"
    pass_when: |
      For planar, radial, cylindrical, and region-distance coordinates, a test in
      tests/coordinate.rs evaluates grad_xi at a representative off-singularity point
      and asserts each component agrees with a 3-point central finite difference of
      xi within the stated tolerance (<= 1e-6 absolute); the suite passes.
      (Implemented as inline `#[cfg(test)] mod tests` in coordinate.rs — keeps
      the internal Coordinate type non-public until a later spec wires it.)
    status: verified
  - id: ac-002
    summary: Inside r_guard, xi and grad_xi are finite with bounded gradient norm
    type: code
    evaluator_hint: "test target: tests/coordinate.rs r_guard cases"
    pass_when: |
      A test in tests/coordinate.rs evaluates a radial coordinate AT its centre and a
      cylindrical coordinate ON its axis and asserts xi and every grad_xi component are
      finite (no NaN/Inf) and the gradient norm is bounded (<= 1.0 + tol); the suite passes.
    status: verified
  - id: ac-003
    summary: Planar/radial/cylindrical xi uses minimum-image delta under an active periodic box
    type: code
    evaluator_hint: "test target: tests/coordinate.rs wrap case"
    pass_when: |
      A test in tests/coordinate.rs places a reference point near a box face, queries a
      point across the periodic boundary, and asserts the returned xi equals the
      hand-computed minimum-image (wrapped) value, not the raw distance; the suite
      passes.
    status: verified
---

# Acceptance criteria

- **ac-001 — analytic vs finite-difference ∇ξ.** The chain-rule consumer downstream
  relies on `grad_xi` being the exact gradient of `xi`. This is the load-bearing
  correctness check for the module: for every coordinate flavour, the hand-derived
  analytic gradient must reproduce a central finite difference of the scalar field at
  off-singularity points. Mirror the FD style of
  `assert_gradient_opposes_violation` (tests/restraint.rs:24) but tighten it to a
  strict componentwise analytic-vs-FD comparison of ∇ξ.

- **ac-002 — finite gradient inside r_guard.** r̂ (radial) and ρ⃗/‖ρ⃗‖ (cylindrical)
  are undefined at the singularity, and an unguarded 1/ξ force diverges. The clamp
  must keep ξ and ∇ξ finite and bounded exactly AT the centre / ON the axis, so gencan
  never sees NaN/Inf.

- **ac-003 — minimum-image wrapping.** The reference plane/centre/axis is fixed in the
  box and slabs may straddle the boundary, so planar/radial/cylindrical deltas must be
  wrapped via the existing `delta_vector` helper (src/cell.rs:108) rather than using
  raw Cartesian differences. Verified against a hand-computed wrapped value.
