---
slug: special-bonds-04-residual
created: 2026-09-04
criteria:
  - id: ac-001
    summary: PackResult.intra is filled at assemble from sys.simbox
    type: code
    pass_when: |
      IntraResidual { scored, exempted } has no Default. PackResult has pub intra.
      Pipeline::assemble fills it via from_targets(..., &tables) where tables is
      vec![from_exclusion_depth(3); n] after positions_in_target_order, using
      &sys.simbox. python/src/result.rs is unchanged.
    status: verified
    last_checked: 2026-09-04
  - id: ac-002
    summary: from_targets takes per-target BondDistanceWeights, not a nailed 3
    type: code
    pass_when: |
      from_targets has signature (targets, positions, simbox, tables) -> Self
      and asserts targets.len() == tables.len(). The from_targets body contains
      no from_exclusion_depth(3). There is no second public constructor that
      omits the table.
    status: verified
    last_checked: 2026-09-04
  - id: ac-003
    summary: Linear hexamer goldens are exempted 1.0 Å and scored 4.0 Å
    type: runtime
    pass_when: |
      intra_residual_linear_hexamer_depth3 asserts exempted == 1.0 and scored == 4.0
      on (i,0,0) i=0..5 with a depth-3 table.
    status: verified
    last_checked: 2026-09-04
  - id: ac-004
    summary: Folded 3-4-5 hexamer reports the 1-6 contact as scored min
    type: runtime
    pass_when: |
      scored == 3.0 (pair 0-5) and exempted == 3.0 on the 3-4-5 coordinates
      in the spec Testing strategy.
    status: verified
    last_checked: 2026-09-04
  - id: ac-005
    summary: Depth-1 vs depth-3 on a 5-bead chain classifies 1-3 differently
    type: runtime
    pass_when: |
      Same 5-atom chain yields scored == 4.0 at depth 3 and scored == 2.0 at
      depth 1 (exempted == 1.0 both). A constructor that nails depth 3 fails
      the depth-1 half.
    status: verified
    last_checked: 2026-09-04
  - id: ac-006
    summary: Empty pair class is +infinity; coincident 0 is a real pair
    type: runtime
    pass_when: |
      Bonded dimer: scored.is_infinite(), exempted == 1.0. Bondless coincident
      atoms: scored == 0.0 with no AABB/NeighborList branch.
    status: verified
    last_checked: 2026-09-04
  - id: ac-007
    summary: A second copy is a different molecule and does not enter intra
    type: runtime
    pass_when: |
      count=2 linear hexamer with copy 1 translated 0.25 Å still yields
      exempted == 1.0 and scored == 4.0.
    status: verified
    last_checked: 2026-09-04
  - id: ac-008
    summary: Scored min uses the simulation-box minimum image
    type: runtime
    pass_when: |
      Linear hexamer at x=0,1,2,3,4,9.5 in a 10 Å PBC-x box yields scored == 0.5.
    status: verified
    last_checked: 2026-09-04
  - id: ac-009
    summary: fdist skip of same-molecule pairs is behaviorally unchanged
    type: code
    pass_when: |
      pair_term still returns None for same ibmol/ibtype before any distance.
      git diff on src/objective.rs contains only rustdoc pointing at PackResult.intra.
    status: verified
    last_checked: 2026-09-04
  - id: ac-010
    summary: Per-copy MIC loop; no NeighborList, grid, or fourth walker
    type: code
    pass_when: |
      src/entry/result.rs classifies each copy with an i<j loop calling
      SimBox::shortest_vector_impl; ripgrep finds no NeighborList, AABB,
      OverlapField, latomfirst, or cutoff argument.
    status: verified
    last_checked: 2026-09-04
  - id: ac-011
    summary: Identity-only only for missing template and zero-edge graphs
    type: code
    pass_when: |
      from_coords dimer scores every i!=j pair. Bond (0,9) yields both
      infinities, not the 1-2 length as scored. rustdoc names NotFound /
      Validation as omitted.
    status: verified
    last_checked: 2026-09-04
  - id: ac-012
    summary: Regression example reproduces the linear-hexamer goldens
    type: runtime
    pass_when: |
      regressions/special-bonds-04-residual.md embeds exempted = 1.0 and
      scored = 4.0 and states assemble supplies from_exclusion_depth(3) at
      the call site, with no third-party runtime.
    status: verified
    last_checked: 2026-09-04
  - id: ac-013
    summary: rustdoc states Å, table argument, Target.special_bonds successor
    type: docs
    pass_when: |
      PackResult.intra rustdoc says 04 assemble broadcasts from_exclusion_depth(3)
      and 05 swaps that list builder for Target.special_bonds. from_targets
      rustdoc names the tables argument and InternalTree::from_frame_with_weights
      as the call-site-3 analogue. It does not name resolved_special_bonds.
      cargo doc --no-deps emits no warning.
    status: verified
    last_checked: 2026-09-04
  - id: ac-014
    summary: Full rust check and test suite pass
    type: runtime
    pass_when: |
      cargo fmt -- --check, clippy -D warnings, cargo test -p molcrafts-molpack
      --lib --tests, and cargo doc --no-deps all succeed.
    status: verified
    last_checked: 2026-09-04
out_of_scope:
  - Target.special_bonds (05)
  - Python mirror (06)
  - NeighborList
  - Changing fdist skip
---

# Acceptance — special-bonds-04-residual

Done means PackResult.intra reports intramolecular mins with the table as a
from_targets argument, assemble broadcasts depth 3 at the call site, and the
walk is a per-copy MIC loop with no NeighborList.
