# Regression — special-bonds-05-target

Public-API pin of the per-target intramolecular skip table
(`Target.special_bonds`), the binary compile gate on
`InternalTree::from_frame`, and the C12 skip-set golden. **hard-coded
golden** default table, skip indices, and error slot; this file is not
compiled. Runtime ownership stays in `src/target.rs` and
`src/grow/internal.rs`.

**Provenance.** Captured 2026-09-04 from the `special-bonds-05-target`
spec implementation (`Target::with_special_bonds`,
`InternalTree::from_frame`, `GrowError::NonBinarySpecialBond` as written
in that spec). No third-party runtime: no live literature code, no
subprocess, no Packmol binary, no live molrs oracle. Domain convention:
Cassandra tail-weight table `[0,0,0,1]` ≡ exclusion depth 3
(special-bonds-01 in molrs); slot 2 of a fractional table is the 1-4
weight.

## Public surface

The table is visible only through these public names:

- `Target.special_bonds` — `molrs::BondDistanceWeights`. Default is
  `from_exclusion_depth(3)` (`[0, 0, 0, 1]`), written once in
  `Target::from_parts`. Re-exported at the crate root next to `Element`,
  not in the prelude. No `SpecialBonds` alias.
- `Target::with_special_bonds` — stores the table and returns `Self`.
  Fractional weights are stored here and refused later.
- `InternalTree::from_frame` — the only public constructor. A non-0/1
  weight is `GrowError::NonBinarySpecialBond { index, weight }` before
  exclusions are compiled.
- `InternalTree::exclusions` — per-atom skip set (sorted, self included).

`tree_from_target` and `binary_violation` are crate-private; they are
not a runtime requirement of this pin. There is no
`PackError::NonBinarySpecialBond` (the entry wraps `PackError::Grow`).

## Hard-coded goldens

| Pin | Literal |
|---|---|
| Default table | `[0, 0, 0, 1]` (`Target::new` / `from_coords`) |
| C12 `exclusions(0)` | `[0, 1, 2, 3]` (`InternalTree::from_frame` + `from_exclusion_depth(3)`; self + 1-2/1-3/1-4; atom 4 is scored) |
| Non-binary refusal | `GrowError::NonBinarySpecialBond { index: 2, weight: 0.5 }` on `[0, 0, 0.5, 1]`; Display names `1-4`, the weight, and `with_atom_radius` |

A public `Target` built from coordinates or a frame carries
`special_bonds.as_slice() == [0.0, 0.0, 0.0, 1.0]`. A 12-atom linear
chain built with depth 3 must report `exclusions(0) == [0, 1, 2, 3]`.
A fractional 1-4 weight is stored on the target and refused at
`InternalTree::from_frame` with `index == 2`.

## Owning tests

Do not copy these. Re-run them to refresh the pin:

- `target_default_special_bonds_is_depth_3` — default table is
  `[0, 0, 0, 1]` from both constructors.
- `internal_exclusions_depth` — C12 `from_frame` +
  `from_exclusion_depth(3)` yields `exclusions(0) == [0, 1, 2, 3]`.
- `internal_tree_rejects_non_binary_special_bond` — slot 2 weight 0.5
  is `NonBinarySpecialBond`; Display contains `1-4`, `0.5`, and
  `with_atom_radius`.

```text
cargo test -p molcrafts-molpack --lib --tests -- target_default_special_bonds_is_depth_3
cargo test -p molcrafts-molpack --lib --tests -- internal_exclusions_depth
cargo test -p molcrafts-molpack --lib --tests -- internal_tree_rejects_non_binary_special_bond
```
