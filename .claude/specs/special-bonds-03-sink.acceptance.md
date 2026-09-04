---
slug: special-bonds-03-sink
created: 2026-09-04
criteria:
  - id: ac-001
    summary: frame_positions is pub(crate) in src/template.rs with a local error
    type: runtime
    pass_when: |
      src/template.rs defines pub(crate) fn frame_positions -> Result<Vec<[F; 3]>,
      FramePositionsError::NoAtomsBlock>. It does not mention MolRsError.
      src/frame.rs does not define frame_positions. src/lib.rs does not re-export
      it. In-module tests round-trip a zigzag frame and refuse missing atoms / missing z.
    status: verified
    last_checked: 2026-09-04
  - id: ac-002
    summary: molpack no longer defines Topology or TopologyError
    type: code
    pass_when: |
      src/topology.rs and tests/topology.rs are absent; src/lib.rs has no pub mod
      topology and no pub use of Topology, TopologyError, or frame_positions;
      there is no pub use molrs::Topology alias.
    status: verified
    last_checked: 2026-09-04
  - id: ac-003
    summary: Bondless templates fail as GrowError::NoBonds
    type: runtime
    pass_when: |
      grow_rejects_template_without_bonds matches GrowError::NoBonds (not Topology
      or MolRs). A 2-atom bondless frame is NoBonds, not TemplateTooSmall or
      Disconnected. Display equals the pre-spec NoBonds sentence.
    status: verified
    last_checked: 2026-09-04
  - id: ac-004
    summary: Disconnected uses n_components; rings cannot build a tree
    type: runtime
    pass_when: |
      from_frame_with_weights returns Disconnected for an isolated atom and for
      two disjoint 3-atom chains. A ring template is RingTemplate and does not
      construct InternalTree. bfs_order no longer returns Disconnected.
    status: verified
    last_checked: 2026-09-04
  - id: ac-005
    summary: Frame-read failures are named GrowError variants, not MolRs(String)
    type: runtime
    pass_when: |
      No atoms / missing z → GrowError::NoAtomsBlock. OOB bond → BondOutOfRange
      { a, b, n } from the frame. GrowError has no MolRs variant and does not
      contain MolRsError.
    status: verified
    last_checked: 2026-09-04
  - id: ac-006
    summary: InternalTree has no from_frame; default-3 lives only on GrowConfig
    type: code
    pass_when: |
      from_frame_with_weights is the only constructor. No from_frame,
      from_frame_with_depth, or DEFAULT_EXCLUDE_BONDS. GrowStage uses
      from_exclusion_depth(config.exclusion_depth). Lattice sites pass
      from_exclusion_depth(3) with a lattice-v1 AA rustdoc note.
      topology_for_growth is the only grow-module caller of crate::template.
      internal.rs and decorate.rs import neither crate::frame nor crate::template.
      BondDistanceWeights comes from molrs core, not molrs::ff.
    status: verified
    last_checked: 2026-09-04
  - id: ac-007
    summary: Default-table growth goldens stay bitwise unedited
    type: runtime
    pass_when: |
      The named tests/grow.rs bitwise pins and tests/pipeline.rs stay green with
      unedited assertion bodies. internal_exclusions_depth still asserts
      exclusions(0) == [0, 1, 2, 3].
    status: verified
    last_checked: 2026-09-04
  - id: ac-008
    summary: Regression example pins public-API goldens
    type: runtime
    pass_when: |
      regressions/special-bonds-03-sink.md names NoBonds Display and C12
      exclusions(0) == [0,1,2,3] under from_frame_with_weights +
      from_exclusion_depth(3), with no third-party runtime.
    status: verified
    last_checked: 2026-09-04
  - id: ac-009
    summary: Architecture map follows the sink
    type: docs
    pass_when: |
      docs/architecture.md has no topology.rs CSR leaf; frame_positions is
      documented on template.rs. dg-refine.md says refine consumes molrs::Topology
      directly, never topology_for_growth. cargo doc --no-deps has 0 warnings
      from these symbols.
    status: pending
  - id: ac-010
    summary: Solver seam still has no molrs::ff dependency
    type: code
    pass_when: |
      src/template.rs and src/grow/** do not import molrs::ff or call
      ff::extract_coords.
    status: verified
    last_checked: 2026-09-04
  - id: ac-011
    summary: Growth-owned neighbor-order and root-inclusive pins live in tests/grow.rs
    type: runtime
    pass_when: |
      tests/grow.rs contains a shuffled-bond fixture on which InternalTree
      decomposes and exclusions(root) includes the root. No test in this crate
      calls molrs::Topology::neighbors as if it were a molpack type.
    status: verified
    last_checked: 2026-09-04
out_of_scope:
  - special-bonds-04 residual
  - special-bonds-05 Target table
  - special-bonds-06 Python
  - special-bonds-02 ladder
  - automating MOLRS_GIT_REF
  - GrowError::MolRs / MolRsError on the packing surface
---

# Acceptance — special-bonds-03-sink

Done means molpack Topology is gone, growth reads molrs::Topology, geometry
lives in src/template.rs, GrowError read failures are matchable, InternalTree
has no default-3 adapter, and bitwise goldens stay green.
