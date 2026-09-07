---
slug: stl-region-01-region
spec: stl-region-01-region
created: 2026-09-05
criteria:
  - id: ac-001
    summary: region.rs moves to region/mod.rs; stl module is private
    type: code
    pass_when: |
      src/region.rs is absent. src/region/mod.rs defines Region and the
      existing *Region types, has private mod stl, and pub-uses StlRegion
      and StlError. Existing combinator tests remain in mod.rs.
    status: pending
  - id: ac-002
    summary: StlRegion impls Region; no StlRestraint type
    type: code
    pass_when: |
      src/region/stl.rs defines pub struct StlRegion and impl Region
      (contains, signed_distance, signed_distance_grad, bounding_box,
      declared_cell=None). No StlRestraint symbol under src/.
    status: pending
  - id: ac-003
    summary: Constructor stack file→bytes→triangles; triangles already Å
    type: code
    pass_when: |
      from_file(path, scale) and from_bytes(bytes, scale) return
      Result<StlRegion, StlError> and apply scale at load.
      from_triangles takes Å triangles and has no scale argument.
      All three run the same watertight/degenerate/empty gates.
    status: pending
  - id: ac-004
    summary: Named rejects empty, scale, degenerate, open mesh
    type: runtime
    pass_when: |
      In-module tests match StlError::Empty, InvalidScale,
      DegenerateTriangle, NotWatertight; a 12-triangle unit cube is Ok.
    status: pending
  - id: ac-005
    summary: d < 1e-9 Å is inside with φ=0
    type: runtime
    pass_when: |
      Unit cube point [1.0, 0.5, 0.5] has contains==true and
      signed_distance.abs() < 1e-9.
    status: pending
  - id: ac-006
    summary: Off-surface membership is half-open MT even-odd
    type: runtime
    pass_when: |
      Nested cavity: outer [-2,2]^3 plus inner [-1,1]^3; origin
      contains==false φ==+1; [1.5,0,0] contains==true φ==-0.5;
      [3,0,0] contains==false φ==+1.
    status: pending
  - id: ac-007
    summary: φ magnitude is Eberly closest-point; negative inside
    type: runtime
    pass_when: |
      Unit cube: φ([0.5,0.5,0.5])==-0.5; φ([2,0.5,0.5])==+1;
      φ([3,3,3])==2*sqrt(3) within 1e-12.
    status: pending
  - id: ac-008
    summary: Mesh storage is private; no public Arc
    type: code
    pass_when: |
      Triangle buffer is a private field. No public signature takes or
      returns Arc of the mesh. Weld policy is bitwise identity after scale.
    status: pending
  - id: ac-009
    summary: StlError is local; PackError unchanged
    type: code
    pass_when: |
      StlError is Clone with Io, Parse, Empty, InvalidScale,
      DegenerateTriangle, NotWatertight. src/error.rs has no STL variant.
    status: pending
  - id: ac-010
    summary: Crate root re-exports StlRegion and StlError
    type: code
    pass_when: |
      src/lib.rs pub-uses both at crate root and in prelude; rustdoc table
      names StlRegion.
    status: pending
  - id: ac-011
    summary: RegionRestraint lift uses existing generic
    type: runtime
    pass_when: |
      RegionRestraint(cube) and cube.into_restraint(): f at centre == 0;
      f at [2,0.5,0.5] == 1.0 with scale2=1.
    status: pending
  - id: ac-012
    summary: ASCII and length-matched binary STL parse
    type: runtime
    pass_when: |
      Hard-coded ASCII and binary (84+12*50, including solid-prefixed
      header still binary by length) match from_triangles φ within 1e-9.
    status: pending
  - id: ac-013
    summary: grow, gencan, python untouched
    type: code
    pass_when: |
      Diff does not modify src/grow/, src/gencan/, or python/.
    status: pending
  - id: ac-014
    summary: Regression pins public cube SDF
    type: runtime
    pass_when: |
      regressions/stl-region-01-region.md uses StlRegion::from_triangles
      on a 12-triangle unit cube and hard-codes centre φ==-0.5 and
      [2,0.5,0.5] φ==+1.0.
    status: pending
out_of_scope:
  - python/ and docs/python/ (02-bind)
  - grow/gencan, BVH, snap-ε weld, winding number, PBC mesh
  - StlRestraint, PackError STL variant
---

# Acceptance — stl-region-01-region

Done = watertight mesh is a `Region` with even-odd + Eberly SDF, lifted only via existing `RegionRestraint`.
