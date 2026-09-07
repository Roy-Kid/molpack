---
slug: lattice-branch-saw-02-docs
spec: lattice-branch-saw-02-docs
created: 2026-09-05
criteria:
  - id: ac-001
    summary: test_lattice_grows_star runs LatticeGrow at 2.0 Å without push-off
    type: runtime
    pass_when: |
      python/tests/test_pack_peo_topo.py::TestTinyPack.test_lattice_grows_star
      calls LatticeGrow.run on make_star(2) with with_tolerance(2.0), does not
      raise, asserts grown.natoms == star.n_atoms, and contains no
      GenCanPack, with_restart, C3, or BFS.
    status: verified
    last_checked: 2026-09-05
  - id: ac-002
    summary: test_lattice_refuses_ring matches ValueError ring
    type: runtime
    pass_when: |
      python/tests/test_pack_peo_topo.py::TestNamedErrors.test_lattice_refuses_ring
      calls LatticeGrow.run on make_ring(...) and pytest.raises(ValueError,
      match="ring").
    status: verified
    last_checked: 2026-09-05
  - id: ac-003
    summary: Keep two CbmcGrow star tests as peer-solver witnesses
    type: code
    pass_when: |
      test_cbmc_grows_one_star and test_auhl_two_stars_push_off_converges
      still exist in python/tests/test_pack_peo_topo.py; the auhl test still
      hard-codes CbmcGrow with_tolerance(0.6) then GenCanPack.with_restart
      at 2.0 Å.
    status: verified
    last_checked: 2026-09-05
  - id: ac-004
    summary: Remove test_lattice_refuses_star
    type: code
    pass_when: |
      python/tests/test_pack_peo_topo.py has no test_lattice_refuses_star
      and no pytest.raises expecting LatticeGrow to reject a star.
    status: verified
    last_checked: 2026-09-05
  - id: ac-005
    summary: pack_star is LatticeGrow@2.0 then with_restart@2.0
    type: code
    pass_when: |
      python/examples/pack_peo_topo.py pack_star builds make_star, Target
      via _target (PEO_H_RADIUS), LatticeGrow with_tolerance(2.0), then
      GenCanPack.with_restart with_tolerance(2.0); it contains no
      try/except around LatticeGrow, no CbmcGrow call, no print of Auhl,
      and no PEO_GROW_TOL / 0.6 on that path.
    status: verified
    last_checked: 2026-09-05
  - id: ac-006
    summary: LatticeGrow API page states trees/RingTemplate/non-tetrahedral
    type: docs
    pass_when: |
      docs/python/api-reference.md LatticeGrow section states that trees
      (including branched) are accepted, a cycle raises named RingTemplate
      ValueError, and non-tetrahedral templates stay named rejections; it
      does not contain v1, linear-only, C3, BFS, or “stars use lattice”.
    status: verified
    last_checked: 2026-09-05
  - id: ac-007
    summary: growth.md walks pack_star as an explicit caller pick
    type: docs
    pass_when: |
      docs/python/guide/growth.md “Branched trees and rings” presents
      pack_peo_topo.pack_star as LatticeGrow then caller-side
      with_restart, keeps CbmcGrow as a peer tree grower, states both
      growers raise RingTemplate on cycles and pack_ring picks rigid
      GenCanPack, and does not say “stars are packed with LatticeGrow”.
    status: verified
    last_checked: 2026-09-05
  - id: ac-008
    summary: Regression markdown pins the public-contract literals
    type: runtime
    pass_when: |
      regressions/lattice-branch-saw-02-docs.md exists and hard-codes:
      example pick LatticeGrow then GenCanPack.with_restart at 2.0 Å
      (not Auhl, not CbmcGrow fallback); trees including branched
      accepted; cycles RingTemplate; forbids StarGrow,
      pack_star_lattice, C3/BFS as public contract, linear-only/v1 on
      LatticeGrow, and packmol in identifiers.
    status: verified
    last_checked: 2026-09-05
  - id: ac-009
    summary: Owned files forbid StarGrow, C3/BFS contract, packmol ids
    type: code
    pass_when: |
      The five owned files contain no StarGrow, no pack_star_lattice,
      no C3/BFS presented as public contract, no packmol identifier,
      and no LatticeGrow “linear-only” / “v1” claim.
    status: verified
    last_checked: 2026-09-05
  - id: ac-010
    summary: PEO_GROW_TOL remains on remaining CbmcGrow auhl paths
    type: code
    pass_when: |
      test_auhl_two_stars_push_off_converges still uses 0.6 Å on
      CbmcGrow; examples/pack_peo/main.rs auhl still reads PEO_GROW_TOL
      default 0.6; pack_star does not.
    status: verified
    last_checked: 2026-09-05
out_of_scope:
  - Rust / PyO3 / 01-walk internals (C3, BFS)
  - New public types (StarGrow, pack_star_lattice)
  - Changing CbmcGrow Auhl paths or pack_ring’s GenCanPack pick
  - Binding stars to LatticeGrow in docs/examples.md
  - Lattice ring closure, non-tetrahedral lattice mapping, Pipeline wrapper
---

# Acceptance — lattice-branch-saw-02-docs

本段「完成」= 星形 Python 示例是调用方对 `LatticeGrow` @ 2.0 Å 再 `with_restart` @ 2.0 Å 的一次 pick；文档公布树/环/非四面体验收集而不把拓扑写成算法；测试钉格相收星、拒环，并留下两条 `CbmcGrow` 星形见证。
