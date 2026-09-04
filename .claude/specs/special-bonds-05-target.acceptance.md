---
slug: special-bonds-05-target
created: 2026-09-04
criteria:
  - id: ac-001
    summary: Default special-bonds table is depth 3 ([0,0,0,1])
    type: code
    pass_when: |
      Target::new / from_coords has special_bonds.as_slice() == [0.0, 0.0, 0.0, 1.0].
      from_parts is the only production site that writes from_exclusion_depth(3).
      ripgrep finds no DEFAULT_EXCLUDE_BONDS, from_frame_with_depth, or
      GrowConfig.exclusion_depth in src/.
    status: verified
    last_checked: 2026-09-04
  - id: ac-002
    summary: with_special_bonds stores BondDistanceWeights and returns Self
    type: code
    pass_when: |
      with_special_bonds(from_exclusion_depth(2)) stores [0.0, 0.0, 1.0] and
      returns Self, not Result. There is no resolved_special_bonds and no
      Target::with_exclusion_depth. Crate root re-exports BondDistanceWeights
      like Element with no SpecialBonds alias.
    status: verified
    last_checked: 2026-09-04
  - id: ac-003
    summary: Fractional weights error only at InternalTree::from_frame
    type: code
    pass_when: |
      with_special_bonds of [0,0,0.5,1] returns Self and stores 0.5.
      InternalTree::from_frame returns GrowError::NonBinarySpecialBond { index: 2,
      weight: 0.5 }; Display contains 1-4, 0.5, and with_atom_radius.
      PackError has no NonBinarySpecialBond variant.
    status: verified
    last_checked: 2026-09-04
  - id: ac-004
    summary: Engine exclusion-depth knob is deleted with no alias
    type: code
    pass_when: |
      ripgrep over src/ finds no with_exclusion_depth and no exclusion_depth field.
      grow_config_builder_chain no longer calls it.
    status: verified
    last_checked: 2026-09-04
  - id: ac-005
    summary: from_frame is the binary gate; tree_from_target is pub(crate)
    type: code
    pass_when: |
      InternalTree::from_frame requires &BondDistanceWeights and calls
      binary_violation. tree_from_target is pub(crate). GrowStage,
      LatticeStage, and LatticeGrow::validate_targets all call it.
    status: verified
    last_checked: 2026-09-04
  - id: ac-006
    summary: Per-target tables are not min-folded; two-table GrowStage test
    type: code
    pass_when: |
      tests/grow.rs has a GrowStage::from_targets fixture with two templates
      and distinct BondDistanceWeights whose exclusions and IntraResidual
      classes differ. KG tests set Target::with_special_bonds(from_exclusion_depth(2))
      and keep existing assertions. assemble reads Target.special_bonds.
    status: verified
    last_checked: 2026-09-04
  - id: ac-007
    summary: C12 exclusion goldens and field skip-set shape stay bitwise
    type: runtime
    pass_when: |
      internal_exclusions_depth still asserts exclusions(0) == [0,1,2,3].
      OverlapField::probe/nearest still take excluded: &[u32]. src/grow/field.rs
      is not in this spec's Files.
    status: verified
    last_checked: 2026-09-04
  - id: ac-008
    summary: Python CbmcGrow.with_exclusion_depth is unhooked in 05
    type: code
    pass_when: |
      python/src/entry.rs has no exclusion_depth field or method.
      molpack.pyi has no CbmcGrow.with_exclusion_depth.
      cargo check --manifest-path python/Cargo.toml succeeds.
      Python Target.with_special_bonds is absent (06).
    status: verified
    last_checked: 2026-09-04
  - id: ac-009
    summary: rustdoc states binary default, H-radius advice, not depth 5
    type: docs
    pass_when: |
      with_special_bonds rustdoc states this is not ForceField::special_bonds,
      default ≡ from_exclusion_depth(3), all-atom advice is with_atom_radius
      on explicit H at depth 3, and [0,0,0,0,0,1] is not recommended.
    status: verified
    last_checked: 2026-09-04
  - id: ac-010
    summary: DRAFT specs delete leftover exclusion-depth knobs
    type: docs
    pass_when: |
      dg-refine.md has no intra_exclusion_depth. dg-refine.acceptance.md ac-003
      uses two Targets with distinct BondDistanceWeights, not a run-level depth.
      grow-axes.md has no WalkGrow::with_exclusion_depth.
    status: verified
    last_checked: 2026-09-04
  - id: ac-011
    summary: Regression example reproduces default table and C12 goldens
    type: runtime
    pass_when: |
      regressions/special-bonds-05-target.md embeds default [0,0,0,1],
      C12 exclusions(0) == [0,1,2,3], and NonBinarySpecialBond at index 2,
      with no third-party runtime.
    status: verified
    last_checked: 2026-09-04
  - id: ac-012
    summary: Full Rust fmt, clippy, fast-tier tests, and Python crate check
    type: runtime
    pass_when: |
      cargo fmt -- --check, clippy -D warnings, cargo test -p molcrafts-molpack
      --lib --tests, and cargo check --manifest-path python/Cargo.toml are green.
    status: verified
    last_checked: 2026-09-04
out_of_scope:
  - Python Target.with_special_bonds and docs (06)
  - PackError::NonBinarySpecialBond
  - field.rs signature change
  - CLI / .inp
---

# Acceptance — special-bonds-05-target

Done means the table lives on Target, growth compiles it at InternalTree::from_frame
with a named GrowError, the engine knob is gone in Rust and Python, and 04's
residual reads the per-target table.
