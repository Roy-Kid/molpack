---
slug: stl-region-02-bind
spec: stl-region-02-bind
created: 2026-09-05
criteria:
  - id: ac-001
    summary: from_file loads a watertight 1 Å cube
    type: runtime
    pass_when: |
      python/tests/test_region.py TestStlRegion.from_file succeeds on a
      tmp_path ASCII watertight unit cube with default scale and scale=1.0.
    status: pending
  - id: ac-002
    summary: from_file named-rejects missing and leaky meshes
    type: runtime
    pass_when: |
      OSError for missing path; ValueError for a non-watertight file;
      TypeError for StlRegion(); ValueError for scale=0.
    status: pending
  - id: ac-003
    summary: with_restraint and with_atom_restraint accept StlRegion
    type: runtime
    pass_when: |
      Target.with_restraint(stl) and with_atom_restraint([0], stl) do not
      TypeError; tests do not call GenCanPack/CbmcGrow/LatticeGrow.run.
    status: pending
  - id: ac-004
    summary: PyStlRegion lives in python/src/region.rs
    type: code
    pass_when: |
      python/src/region.rs defines pyclass StlRegion wrapping
      molpack::StlRegion; constraint.rs has no STL parser.
    status: pending
  - id: ac-005
    summary: Extractors lift via RegionRestraint, not collective
    type: code
    pass_when: |
      try_atom_builtin and extract_restraint wrap RegionRestraint(inner);
      extract_collective_restraint has no StlRegion arm.
    status: pending
  - id: ac-006
    summary: Wheel exports StlRegion; no StlRestraint
    type: code
    pass_when: |
      add_class, __all__, and molpack.pyi expose StlRegion.from_file;
      no StlRestraint or packmol identifier in owned files.
    status: pending
  - id: ac-007
    summary: Docs name the three-entry split, not a Force-only story
    type: docs
    pass_when: |
      restraints.md and api-reference.md state Å and centre-only, and
      that GenCanPack is soft, CbmcGrow propose is hard, LatticeGrow
      scores frest after decoration; they do not infer the algorithm
      from the mesh.
    status: pending
  - id: ac-008
    summary: Regression pins the public bind contract
    type: runtime
    pass_when: |
      regressions/stl-region-02-bind.md hard-codes from_file scale=1.0 Å,
      per-atom RegionRestraint lift, three-entry split, and forbids
      StlRestraint.
    status: pending
out_of_scope:
  - src/region/stl.rs, grow, gencan
  - StlRestraint, collective routing, .inp keywords
---

# Acceptance — stl-region-02-bind

Done = Python `StlRegion.from_file` attaches per-atom via `RegionRestraint`; docs tell the three-entry truth.
