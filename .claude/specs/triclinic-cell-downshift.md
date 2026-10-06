---
title: triclinic-cell-downshift — pack into arbitrary lattices by routing all cell/PBC geometry through molrs
status: draft
created: 2026-07-25
---

# triclinic-cell-downshift — pack into arbitrary lattices

## Summary

Delete `molpack/src/cell.rs` (a direct port of Packmol's `cell_indexing.f90` +
`pbc.f90`, orthorhombic and axis-aligned by construction) and route every
cell-assignment, stencil and minimum-image operation through molrs
`SimBox` + `CellGrid`. Packing then works in **any** lattice — hexagonal,
monoclinic, triclinic — with per-axis periodicity, which Packmol cannot do at
all.

The linked-cell *data structures* in `PackContext` (`latomfirst`, `latomnext`,
`lcellfirst`, `neighbor_cells_f/g`, …) and the objective's pair-loop layout are
**kept unchanged**. Only the geometry underneath them is replaced. This keeps
the hot loop byte-for-byte the same shape, so the orthorhombic path stays a
provable no-op and the performance risk is confined to cell assignment.

## Domain basis

- Packmol supports **orthorhombic PBC only** (user guide, v20.15.0: *"Packmol
  supports Orthorhombic periodic boundary conditions"*). Non-orthorhombic cells
  are unsupported; the maintainer confirms in m3g/packmol#120 that it is *"still
  not possible"* and defers it to Packmol.jl, which at the time of writing is
  still only a runner for the Fortran binary.
- Consumers work around this by packing into an orthorhombic bounding box and
  cropping, which loses density near the boundary (SCM/AMS documents that for
  non-orthorhombic cells the result is *approximate* and the density is
  *typically lower than requested*), or by building a rectangular supercell,
  which inflates the atom count (a hexagonal cell needs a √3 rectangularisation,
  ≥ 2× the atoms) and therefore the downstream MD cost.
- m3g/packmol#121 records why this is architectural rather than incidental: in
  Packmol every constraint is evaluated in the reference frame of the "first
  box", which is also the bounding box of the cell-list partitioning, so any
  atom outside it is dropped from the minimum-distance computation. The
  maintainer's own open question — *"what does 'above plane' mean in the
  presence of PBCs?"* — is answered explicitly by this spec (see Design).
- Once cells are indexed in fractional space the stencil in index space is the
  same for orthorhombic and triclinic lattices; the lattice re-enters only
  through the minimum-image displacement. Triclinic support therefore requires
  correct cell sizing (`nearest_plane_distance`) and correct stencils, both of
  which land in molrs under spec `cell-grid-api`.

## Design

### Cell geometry state

`PackContext` today carries `ncells`, `cell_length`, `pbc_length`, `pbc_min`,
`pbc_periodic`. These are replaced by a single `SimBox` plus a `CellGrid`:

```rust
pub struct PackContext {
    pub simbox: SimBox,        // lattice matrix + origin + per-axis pbc
    pub grid: CellGrid,        // celldim + pbc, from CellGrid::for_cutoff
    // linked-cell lists unchanged: latomfirst / latomnext / lcellfirst / ...
}
```

Call-site mapping:

| today | after |
|---|---|
| `cell::setcell(pos, pbc_min, pbc_length, cell_length, ncells, periodic)` | `grid.cell3(&simbox, pos)` |
| `cell::index_cell(cell, ncells)` | `grid.flat(cell)` |
| `cell::icell_to_cell(icell, ncells)` | `grid.unflat(icell)` |
| `cell::cell_ind(idx, ncells)` | folded into `CellGrid` stencils |
| `cell::delta_vector(xi, xj, pbc_length, periodic)` | `simbox.shortest_vector_impl(xi, xj)` |

`neighbor_cells_f` / `neighbor_cells_g` (the precomputed 13 forward neighbours
per cell) are rebuilt from `CellGrid::stencil_forward`, which also fixes the
known half-stencil miscount when a periodic axis holds fewer than 3 cells.

### Constraint semantics under PBC

This is the answer to m3g/packmol#121, and it is a design decision, not an
implementation detail:

- **Region restraints are expressed in fractional coordinates.** `InsideBox`
  gains a fractional sibling `InsideCell` meaning "inside the primitive cell";
  the boolean region algebra (`And`/`Or`/`Not`) composes unchanged.
- **A half-space restraint on a periodic axis is rejected, not reinterpreted.**
  `AbovePlane` / `BelowPlane` whose normal has a non-zero component along a
  periodic lattice direction returns `PackError` at build time with a message
  naming the offending axis. A plane has no meaning on a periodic axis; silently
  evaluating it in the "first box" frame — Packmol's behaviour — produces a
  constraint whose satisfaction depends on where the user happened to put the
  origin.
- **Fixed molecules may cross the periodic boundary.** Because the cell
  assignment clamps on non-periodic axes and wraps on periodic ones, a fixed
  slab or protein straddling a periodic face is still indexed correctly. Packmol
  forbids this (`initial.f90` hard error).

### Builder / script surface

```rust
Molpack::with_cell(lengths: [F; 3], angles_deg: [F; 3], pbc: [bool; 3])
Molpack::with_cell_matrix(h: [[F; 3]; 3], origin: [F; 3], pbc: [bool; 3])
Molpack::with_periodic_box(min, max)   // retained: orthorhombic shorthand
```

`with_periodic_box` currently forces all three axes periodic, which makes a slab
pack (xy periodic, z confined) stall badly; it takes the per-axis `pbc` from the
new API instead of hard-coding `[true; 3]`. The `.inp` script grows a `cell`
keyword (`cell a b c alpha beta gamma` and `cell_matrix …`) alongside the
existing `pbc`.

### What is deliberately not changed

The pair loop keeps Packmol's linked-list traversal and its 13-forward-neighbour
precompute. Swapping in molrs `LinkCell`'s counting-sorted layout and
`visit_pairs` is a separate, later question: it changes pair emission order and
therefore floating-point summation order, so it cannot be validated by the
"orthorhombic path is unchanged" argument this spec relies on.

## Files

- `molpack/src/cell.rs` — **deleted**
- `molpack/src/context/pack_context.rs` — `SimBox` + `CellGrid` fields; rebuild
  `neighbor_cells_f/g` from `CellGrid::stencil_forward`
- `molpack/src/objective.rs` — cell assignment and minimum image via molrs
- `molpack/src/initial.rs` — random placement drawn in fractional coordinates;
  `avoid_overlap` stencil scan via `CellGrid`
- `molpack/src/restraint/mod.rs`, `molpack/src/restraint/geometric/bounded.rs` —
  fractional `InsideCell`; periodic-axis rejection for plane restraints
- `molpack/src/region.rs` — `CellRestraint` in the `And`/`Or`/`Not` algebra
- `molpack/src/packer.rs`, `molpack/src/handler.rs`, `molpack/src/script/build.rs`
  — builder + script plumbing
- `molpack/src/context/work_buffers.rs` — buffer sizing from `CellGrid::n_cells`
- ~~`molpack/benches/pack_end_to_end.rs`, `molpack/benches/pair_kernel.rs` —
  triclinic variants added~~ (`benches/` deleted 2026-09-20; no replacement yet)
- ~~`molpack/tests/` — new `triclinic.rs`; existing `packer.rs`,
  `examples_batch.rs`, `gradient.rs` must pass unchanged~~ (`tests/` deleted
  2026-09-20; coverage is in-module, e.g. `entry::setup` periodic-declaration tests)

## Tasks

1. Carry `SimBox` + `CellGrid` in `PackContext`; mechanical call-site swap;
   orthorhombic behaviour byte-identical.
2. Rebuild `neighbor_cells_f/g` from `CellGrid::stencil_forward`; delete
   `src/cell.rs`.
3. Fractional `InsideCell` region + boolean-algebra integration.
4. Plane-restraint periodic-axis rejection with a named-axis error.
5. Builder + `.inp` `cell` keyword; fix `with_periodic_box` per-axis pbc.
6. Triclinic tests + triclinic bench variants.

## Testing

**Orthorhombic no-op.** The five shipped Packmol examples (`mixture`,
`interface`, `bilayer`, `spherical`, `solvprotein`) run at a fixed seed and
produce the same final objective and the same validation report as before the
change. This is the Case 0 regression of the paper and the guard for the whole
refactor.

**Triclinic correctness.** For a hexagonal cell (a = b, γ = 120°) and a
strongly tilted triclinic cell, packing a simple solvent to a stated tolerance
must satisfy that tolerance under the **true** minimum image. The oracle is an
independent 27-image brute-force scan (not the packer's own cell list): no
inter-molecular atom pair may be closer than the tolerance.

**Density.** Packing N molecules into a hexagonal cell of volume V yields
density N/V within the packer's own precision, and is compared against the
orthorhombic-bounding-box workaround (pack + crop) on the same lattice, which is
expected to fall short — this comparison is the paper figure.

**Mixed periodicity.** A slab (xy periodic, z confined by `InsideCell` on the
z fraction only) converges without the stall currently caused by
`with_periodic_box` forcing all axes periodic.

## Out of scope

- Replacing the linked-list pair traversal with molrs `LinkCell::visit_pairs`.
- Force-field-based relaxation; nothing in this spec requires the `ff` feature.
- Collective/distribution restraints — separate spec
  `collective-com-restraints`, which depends on this one only for the
  reciprocal-lattice vectors of a general cell.
- Any performance claim. Timing appears in the paper as a compatibility note
  only; the gate here is "no catastrophic regression", not "faster".
