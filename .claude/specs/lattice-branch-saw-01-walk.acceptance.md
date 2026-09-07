---
slug: lattice-branch-saw-01-walk
spec: lattice-branch-saw-01-walk
created: 2026-09-05
criteria:
  - id: ac-001
    summary: Backbone owns parent/children/align; no WalkTree
    type: code
    pass_when: |
      src/grow/lattice/decorate.rs defines pub(crate) Backbone with atoms,
      parent: Vec<Option<usize>>, children: Vec<Vec<usize>>, align: [usize; 3],
      parent_bond, and follows: Vec<Option<(usize, F)>> parallel to atoms;
      no WalkTree or BranchedBackbone type exists in src/grow/lattice/.
    status: verified
    last_checked: 2026-09-05
  - id: ac-002
    summary: grow_walk borrows parent/children/follows; no decorate import
    type: code
    pass_when: |
      grow_walk and forced_zigzag in src/grow/lattice/saw.rs take
      parent: &[Option<usize>], children: &[Vec<usize>], and
      follows: &[Option<(usize, F)>] (or equivalent); they do not import
      decorate::Backbone and have no hydrogen-mask argument.
    status: verified
    last_checked: 2026-09-05
  - id: ac-003
    summary: Star-decision stack replaces tried Vec; one free draw per var
    type: code
    pass_when: |
      src/grow/lattice/saw.rs contains no tried: Vec<Vec<…>> local. At a
      parent, only the hooked child of an InternalTree variable is
      Rosenbluth-sampled; remaining slots (same-v followers and
      follows=None d=1) occupy the tetrahedral neighbour nearest
      follows.offset (not a second draw). Recoil of that decision
      releases those sites. The 5-atom centre+4-leaf public star has
      n_vars=0 and does not pin Rosenbluth-once.
    status: verified
    last_checked: 2026-09-05
  - id: ac-004
    summary: k=1 linear chains use the same grow_walk
    type: runtime
    pass_when: |
      cargo test -p molcrafts-molpack --lib --tests --
      lattice_grow_bead_chain_constructive is green; there is no separate
      linear-only walk function in src/grow/lattice/saw.rs.
    status: verified
    last_checked: 2026-09-05
  - id: ac-005
    summary: parent is InternalTree projected onto not-H atoms
    type: code
    pass_when: |
      analyze_backbone roots at the first not-H atom in InternalTree
      order (seed then step_sites); parent[j] is the nearest heavy on
      that atom's InternalTree parent chain. decorate.rs does not run a
      second heavy-graph diameter BFS to choose the root.
    status: verified
    last_checked: 2026-09-05
  - id: ac-006
    summary: One hook per var uses site offset; offset 0 not required
    type: runtime
    pass_when: |
      analyze_backbone records at most one hook per InternalTree variable
      from a backbone (heavy) Site using that site's follows offset. A
      C/H butane in-module test hooks the heavy when the offset-0 site
      is H. Alignment atoms have no hooks.
    status: verified
    last_checked: 2026-09-05
  - id: ac-007
    summary: sites[j] is a diamond neighbour of sites[parent[j]]
    type: runtime
    pass_when: |
      In-module grow_walk and decorate_chain tests assert that for j not
      the root, Walk.sites[j] is a diamond neighbour of
      Walk.sites[parent[j]], and w[j] is rebuilt from w[parent[j]] along
      parent_bond[j]. No windows(2) / last-site extension remains.
    status: verified
    last_checked: 2026-09-05
  - id: ac-008
    summary: Target dihedral is four points; refs[2] is never a Walk site
    type: code
    pass_when: |
      decorate_chain solves vars from a four-point dihedral: parent-chain
      (ggp, gp, p, j) when parent³ exists, else trial-at-0 on site.refs
      using coords[refs[2]] from place_seed/prefix, never Walk.sites
      indexed by refs[2].
    status: verified
    last_checked: 2026-09-05
  - id: ac-009
    summary: d==0 and d>4 named; d 1..=4 legal; no staged wording
    type: runtime
    pass_when: |
      In-module d==0 rejects as NonTetrahedralTemplate without dropping
      the atom from the walk. tests/grow.rs lattice_grow_rejects_degree_gt_4
      is Err whose message does not contain "branched staged".
      lattice_grow_tetrahedral_star_completes and
      lattice_grow_tetrahedral_comb_completes return Ok.
    status: verified
    last_checked: 2026-09-05
  - id: ac-010
    summary: not-H lives only in analyze_backbone; no Target H API
    type: code
    pass_when: |
      element==H (AA default) and missing element → "X" are documented
      and implemented only in analyze_backbone; saw.rs has no copy of
      that mask; Target gains no hydrogen method.
    status: verified
    last_checked: 2026-09-05
  - id: ac-011
    summary: In-module tests own pub(crate) walk and backbone
    type: runtime
    pass_when: |
      src/grow/lattice/decorate.rs and saw.rs each have #[cfg(test)] mod
      tests covering analyze_backbone/decorate_chain and
      grow_walk/forced_zigzag; cargo test -p molcrafts-molpack --lib
      --tests -- grow::lattice is green. tests/grow.rs does not call
      grow_walk.
    status: verified
    last_checked: 2026-09-05
  - id: ac-012
    summary: LatticeGrow::run completes on tetrahedral star and comb
    type: runtime
    pass_when: |
      tests/grow.rs lattice_grow_tetrahedral_star_completes yields Ok with
      natoms()==10 and lattice_grow_tetrahedral_comb_completes yields Ok
      with natoms()==24; lattice_grow_rejects_branched is absent.
    status: verified
    last_checked: 2026-09-05
  - id: ac-013
    summary: LatticeStage still requires None and guarantees All
    type: runtime
    pass_when: |
      cargo test -p molcrafts-molpack --lib --tests --
      lattice_stage_requires_none_guarantees_all is green.
    status: verified
    last_checked: 2026-09-05
  - id: ac-014
    summary: Crate rustdoc treats tetrahedral trees as first-class
    type: docs
    pass_when: |
      rustdoc on src/grow/lattice/mod.rs, entry.rs LatticeGrow,
      decorate.rs, saw.rs, and GrowError::NonTetrahedralTemplate in
      src/grow/config.rs describe tetrahedral trees of degree <= 4;
      python/ and docs/python/ are untouched; internal.rs is untouched.
    status: verified
    last_checked: 2026-09-05
  - id: ac-015
    summary: Regression pin runs public LatticeGrow on a 5-atom star
    type: runtime
    pass_when: |
      regressions/lattice-branch-saw-01-walk.md shows LatticeGrow::run on
      a 5-atom tetrahedral star, 2 copies, seed 7, box 20 Å, and
      hard-codes Ok, natoms==10, and center–leaf bonds 1.53±1e-6; no
      third-party runtime.
    status: verified
    last_checked: 2026-09-05
  - id: ac-016
    summary: Branch continuations are unused T_STEPS with T_i·T_j == -1
    type: runtime
    pass_when: |
      An in-module grow_walk star test asserts each placed child
      direction is a T_STEPS index unused by the parent and pairwise
      T_i·T_j == -1 (hard-coded lattice identity, no third-party).
    status: verified
    last_checked: 2026-09-05
out_of_scope:
  - python/ and docs/python/ (lattice-branch-saw-02-docs)
  - Target decoration-set API
  - rings, Space trait, lattice MC repair, editing internal.rs
---

# Acceptance — lattice-branch-saw-01-walk

本段「完成」= 四面体重原子树经现有 `LatticeGrow` 在金刚石格上走完（线性为 \(d=2,k=1\)），`parent` 是 InternalTree 投影，每个变量只抽一次方向，装饰挂钩跟真实 `Site::follows` offset，对齐三点不挂钩。Python 站点仍欠 02。
