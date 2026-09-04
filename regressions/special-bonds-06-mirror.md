# Regression — special-bonds-06-mirror

Public-API pin of the Python wheel after 05: the intramolecular skip
table is `Target.with_special_bonds` / `Target.special_bonds`, the
engine knob `CbmcGrow.with_exclusion_depth` is absent, and
`State.intra` forwards `{ scored, exempted }`. **hard-coded
golden** default table and the 2026-09-04 H-radius measurement; this
file is not compiled. Runtime ownership stays in `python/src/target.rs`
and `python/src/result.rs`.

**Provenance.** Captured 2026-09-04 from the `special-bonds-06-mirror`
spec implementation (`PyTarget.with_special_bonds`, `PyIntraResidual`,
and `State.intra` as written in that spec). No third-party
runtime: no live literature code, no subprocess, no Packmol binary, no
live molrs oracle. Domain convention: Cassandra tail-weight table
`[0,0,0,1]` ≡ exclusion depth 3 (special-bonds-01 in molrs). The
2026-09-04 dp5 PEO melt (depth 3, H = 0.85 Å → 142 rounds / 0.3 s) is
commented below and is not re-run from this pin.

## Public surface

The table and residual are visible only through these public names:

- `Target.with_special_bonds(table)` — `list[float]` marshalled by
  `BondDistanceWeights::new` (non-empty, finite, `0 ≤ w ≤ 1`).
  Fractional 0.5 is stored. No `BondDistanceWeights` pyclass.
- `Target.special_bonds` — `inner.special_bonds.as_slice()`. Default
  `[0, 0, 0, 1]`.
- `State.intra` — nested `IntraResidual { scored, exempted }` in
  Å. Forwards the assembled residual; does not recompute from
  positions. No `min_intra_*`.
- `CbmcGrow.with_exclusion_depth` — **absent** (unhooked in 05).

`binary_violation` is not copied into the Python builder; a fractional
table is refused at `CbmcGrow.run` as `PackError::Grow`
(`GrowError::NonBinarySpecialBond`).

## Hard-coded goldens

| Pin | Literal |
|---|---|
| Default table | `[0, 0, 0, 1]` (`Target` constructor; depth 3) |
| Engine knob | `hasattr(CbmcGrow(...), "with_exclusion_depth") is False` |
| IntraResidual | `State.intra.scored` / `.exempted` (Å; empty class `+∞`) |
| H-radius (2026-09-04, dp5 PEO) | depth 3 + H = 0.85 Å → 142 rounds / 0.3 s |

A public `Target` built from a frame carries
`special_bonds == [0.0, 0.0, 0.0, 1.0]`. A fractional 1-4 weight is
stored on the target and refused at `CbmcGrow.run`. All-atom explicit
hydrogen keeps the default table and shrinks H via
`Target.with_atom_radius`; the 2026-09-04 PEO melt is the measured
recipe, not a live assertion.

## Owning tests

Do not copy these. Re-run them to refresh the pin:

- `python/tests/test_target.py` — default table, immutability, CG
  round-trip, empty/non-finite/out-of-range `ValueError`, fractional
  0.5 stores (no `run()`).
- `python/tests/test_pack_result.py` — nested `IntraResidual`
  forwarding; bonded diatomic empty scored class `+∞`.
- `python/tests/test_grow.py` — `hasattr(..., "with_exclusion_depth")
  is False`; fractional 0.5 raises at `CbmcGrow.run`.

```text
uv run --directory python --group dev tox -e py
```
