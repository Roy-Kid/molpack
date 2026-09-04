# Regression — special-bonds-03-sink

Public-API pin of the named bond-less refusal after the molpack Topology
sink, and of the default-3 intramolecular skip set on a linear C12.
**hard-coded golden** Display text and skip indices; this file is not
compiled. Runtime ownership stays in `tests/grow.rs`.

**Provenance.** Captured 2026-09-04 from the `special-bonds-03-sink`
spec implementation (`GrowError::NoBonds` Display and
`InternalTree::from_frame_with_weights` skip sets as written in that
spec). No third-party runtime: no live literature code, no subprocess,
no Packmol binary, no live molrs oracle. Domain convention cited from
the Cassandra tail-weight table `[0,0,0,1]` ≡ exclusion depth 3
(special-bonds-01 in molrs) — recorded here as the literal skip set
`[0, 1, 2, 3]` for atom 0 of a linear C12.

## Public surface

The sink is visible only through these public names:

- `GrowError::NoBonds` — named refusal when the template frame carries
  no bonds (missing or empty graph). `Display` is the sentence below;
  the variant is matchable and does not wrap a molrs error.
- `CbmcGrow::run` — growth entry. A bond-less template fails as
  `PackError::Grow { source: GrowError::NoBonds, .. }`; growth does not
  silently degrade to `GenCanPack`.
- `InternalTree::from_frame_with_weights` — the only public constructor.
  The skip table is `molrs::BondDistanceWeights::from_exclusion_depth`.
- `InternalTree::exclusions` — per-atom skip set (sorted, self included).

`topology_for_growth` and `frame_positions` are crate-private; they are
not a runtime requirement of this pin. `molrs::Topology` is not
re-exported as `molpack::Topology`.

## Hard-coded goldens

| Pin | Literal |
|---|---|
| `GrowError::NoBonds` Display | `the template frame carries no bonds; growth needs the bond graph — pack this target with GenCanPack or supply connectivity` |
| C12 `exclusions(0)` | `[0, 1, 2, 3]` (`InternalTree::from_frame_with_weights` + `BondDistanceWeights::from_exclusion_depth(3)`; self + 1-2/1-3/1-4; atom 4 is scored) |

A public growth run of a template with an atoms block and no bonds block
must fail as `GrowError::NoBonds` with that Display. A 12-atom linear
chain built with depth 3 must report `exclusions(0) == [0, 1, 2, 3]`.

## Owning tests

Do not copy these. Re-run them to refresh the pin:

- `grow_rejects_template_without_bonds` — `CbmcGrow` on a bond-less
  template is `GrowError::NoBonds` with the Display sentence above.
- `internal_exclusions_depth` — C12 `from_frame_with_weights` +
  `from_exclusion_depth(3)` yields `exclusions(0) == [0, 1, 2, 3]`.

```text
cargo test -p molcrafts-molpack --lib --tests -- grow_rejects_template_without_bonds
cargo test -p molcrafts-molpack --lib --tests -- internal_exclusions_depth
```
