# Regression — special-bonds-04-residual

Public-API pin of the intramolecular residual on every engine run:
same-copy **scored** vs **exempted** minima (Å, minimum image).
**hard-coded golden** `exempted = 1.0` and `scored = 4.0` on a linear
hexamer; this file is not compiled. Runtime ownership stays in
`src/entry/result.rs`.

**Provenance.** Captured 2026-09-04 from the `special-bonds-04-residual`
spec implementation (`IntraResidual`, `PackResult.intra`, and
`IntraResidual::from_targets` as written in that spec). Assemble
supplies `from_exclusion_depth(3)` at the call site
(`vec![BondDistanceWeights::from_exclusion_depth(3); n]` in
`Pipeline::assemble`). No third-party runtime: no live literature code,
no subprocess, no Packmol binary, no live molrs oracle. Domain
convention: Packmol's `fdist` is intermolecular only (Martínez et al.,
*J. Comput. Chem.* **30**, 2157 (2009), doi:10.1002/jcc.21224); the
hand-computed hexamer mins below close that silent intramolecular gap.
Cassandra tail-weight table `[0,0,0,1]` ≡ exclusion depth 3
(special-bonds-01 in molrs).

## Public surface

The residual is visible only through these public names:

- `IntraResidual` — `{ scored, exempted }` in Å (minimum image). An
  empty class is `+∞`. There is no `Default`. Re-exported at the crate
  root, not in the prelude.
- `PackResult.intra` — filled at assemble. `fdist` still skips
  same-molecule pairs.
- `IntraResidual::from_targets` — the only public constructor. The skip
  table is the `tables: &[BondDistanceWeights]` argument (one per
  target; `assert_eq!(targets.len(), tables.len())` in release as well).
  There is no table-less overload. Analogous table-passing:
  `InternalTree::from_frame_with_weights`.

Assemble supplies `from_exclusion_depth(3)` at the call site; the
constructor does not nail 3. Spec 05 swaps that list builder for
`Target.special_bonds`; this pin does not name that successor as a
runtime requirement.

The MIC `i < j` loop and `SimBox::shortest_vector_impl` are
crate-private; they are not a runtime requirement of this pin. No
NeighborList, AABB, or OverlapField walk.

## Hard-coded goldens

| Pin | Literal |
|---|---|
| Linear hexamer exempted | `exempted = 1.0` (`(i, 0, 0)` Å, i = 0..5; depth-3 skip of 1-2/1-3/1-4; min is the 1.0 Å bond) |
| Linear hexamer scored | `scored = 4.0` (first scored pair is 1-5 at graph distance 4) |
| Assemble table | `from_exclusion_depth(3)` at the assemble call site (not inside `from_targets`) |

A linear hexamer at `(i, 0, 0)` Å, i = 0..5, classified with
`from_exclusion_depth(3)`, must report `exempted = 1.0` and
`scored = 4.0`. A public engine run fills `PackResult.intra` from that
same depth-3 table because assemble supplies `from_exclusion_depth(3)`
at the call site.

## Owning tests

Do not copy these. Re-run them to refresh the pin:

- `intra_residual_linear_hexamer_depth3` — `IntraResidual::from_targets`
  on a bonded linear hexamer with `from_exclusion_depth(3)` yields
  `exempted == 1.0` and `scored == 4.0`.

```text
cargo test -p molcrafts-molpack --lib --tests -- intra_residual_linear_hexamer_depth3
```
