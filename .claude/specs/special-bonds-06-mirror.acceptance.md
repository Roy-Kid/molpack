---
slug: special-bonds-06-mirror
created: 2026-09-04
criteria:
  - id: ac-001
    summary: Target.with_special_bonds forwards the 05 public field
    type: code
    pass_when: |
      PyTarget.with_special_bonds takes Vec<NpF>; special_bonds getter returns
      inner.special_bonds.as_slice(); default == [0.0, 0.0, 0.0, 1.0];
      BondDistanceWeights is not a Python type; getter does not call
      resolved_special_bonds.
    status: verified
    last_checked: 2026-09-04
  - id: ac-002
    summary: Python builder marshalls BondDistanceWeights::new only
    type: runtime
    pass_when: |
      test_target.py raises ValueError at with_special_bonds for empty,
      non-finite, and w=1.5 tables, stores fractional w=0.5 and round-trips
      it on the getter, and does not call run(). test_grow.py raises at
      CbmcGrow.run on a 0.5 table as PackError.Grow /
      GrowError.NonBinarySpecialBond (1-4 + with_atom_radius). helpers.rs
      is unchanged. python/src/target.rs does not contain binary_violation.
    status: verified
    last_checked: 2026-09-04
  - id: ac-003
    summary: 05 unhook is verified; 06 adds Target.with_special_bonds only
    type: code
    pass_when: |
      This spec's Files list does not include python/src/entry.rs.
      hasattr(CbmcGrow(...), "with_exclusion_depth") is False.
      grep for `fn with_exclusion_depth` / `def with_exclusion_depth` over
      python/src/ and molpack.pyi returns zero hits (the hasattr pin in
      test_grow.py may still mention the name).
      python/src/target.rs defines with_special_bonds.
    status: verified
    last_checked: 2026-09-04
  - id: ac-004
    summary: PackResult.intra mirrors 04 IntraResidual scored/exempted
    type: runtime
    pass_when: |
      PyIntraResidual has scored and exempted; PackResult.intra forwards
      without recomputing; IntraResidual is registered and exported;
      min_intra_scored / min_intra_exempt attributes are absent. A bonded
      diatomic (bonds block present) has infinite scored; _make_tiny_pack
      is used only for names/types.
    status: verified
    last_checked: 2026-09-04
  - id: ac-005
    summary: molpack.pyi lists every Target builder plus intra types
    type: code
    pass_when: |
      molpack.pyi declares with_special_bonds, special_bonds, with_radius,
      with_atom_radius, fscale and short-radius families, IntraResidual,
      PackResult.intra, and does not declare CbmcGrow.with_exclusion_depth.
    status: verified
    last_checked: 2026-09-04
  - id: ac-006
    summary: growth.md AA/CG is Target.with_special_bonds; H radius first-class
    type: docs
    pass_when: |
      growth.md no longer describes exclusion depth as a CbmcGrow knob and
      no longer requires two runs for mixed AA/CG. It tells all-atom callers
      to use with_atom_radius (~0.85 Å) at the default depth-3 table and
      embeds the 2026-09-04 measured comparison. [0,0,0,0,0,1] appears only
      as an escape valve. grep with_exclusion_depth over docs/python/ is 0.
    status: verified
    last_checked: 2026-09-04
  - id: ac-007
    summary: Regression example embeds the public-API goldens
    type: runtime
    pass_when: |
      regressions/special-bonds-06-mirror.md names default [0,0,0,1], no
      with_exclusion_depth, PackResult.intra.scored/.exempted, comments the
      2026-09-04 H-radius table, and imports no third-party packer.
    status: verified
    last_checked: 2026-09-04
  - id: ac-008
    summary: Python gate passes without touching the ABI line
    type: runtime
    pass_when: |
      tox -c python -e py is green. cargo fmt/check on python/Cargo.toml
      succeed. python/src/interop.rs and helpers.rs are unchanged.
    status: verified
    last_checked: 2026-09-04
out_of_scope:
  - Rust Target / GrowConfig API (05)
  - Re-deleting python/src/entry.rs engine knob
  - Recomputing intra or min_intra_* aliases
  - molrs ABI capsules
---

# Acceptance — special-bonds-06-mirror

Done means the wheel tells the same story as Rust after 05 and 04: the table
lives on Target, the engine knob is gone because 05 deleted it, PackResult.intra
reports scored/exempted under those names, Python refuses a fractional table at
the builder, and all-atom PEO is fixed by shrinking explicit-hydrogen radii.
