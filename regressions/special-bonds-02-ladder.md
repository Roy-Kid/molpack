# Regression — special-bonds-02-ladder

Public-API pin of the per-chain softening ladder and the self- vs
inter-chain hard-core split. **hard-coded golden** numbers; this file is
not compiled. Runtime ownership stays in `tests/grow.rs`.

**Provenance.** Captured 2026-09-04 from the `special-bonds-02-ladder`
spec implementation (rung arithmetic and `BlockKind` variants as written
in that spec). No third-party runtime: no live literature code, no
subprocess, no Packmol binary. Domain floor cited from Auhl et al.,
*J. Chem. Phys.* **119**, 12718 (2003), doi:10.1063/1.1628670 — the 0.8σ
push-off bound, recorded here as the literal `0.8`.

## Public surface

The ladder is visible only through these public names:

- `CbmcGrow` — growth entry. `CbmcGrow::with_soften_after(n)` sets the
  dead-end watermark that earns one softening rung.
- `StepInfo.radscale` — on a growth round this **is** the driver's
  hard-core scale. A `Handler::on_step` recorder is the public window;
  it starts at `1.0` and steps down.
- `OverlapField::block_kind` — post-failure query on
  `molpack::grow::field::OverlapField`. Classifies a hard-core hit as
  `BlockKind::SelfBlocked` or `BlockKind::InterChain`. `probe` stays a
  first-hit unit `Blocked` and does not call this.

Private `Chain` clocks are not a runtime requirement of this pin. In
prose they are three readers: `deadend_streak` (cleared on a successful
commit; retract + floor force-place), `deadends_total` (never cleared on
commit; softening watermark with `rungs_earned`), `rungs_earned` (rungs
this chain has taken this run).

## Hard-coded goldens

| Pin | Literal |
|---|---|
| Rung factor | `0.97` (one 3 % shrink per earned rung; at most one per round) |
| Auhl floor | `0.8` (`GrowConfig::new` default `min_hard_scale`; `StepInfo.radscale` never drops below it) |
| Same-molecule core hit | `SelfBlocked` (`BlockKind::SelfBlocked`) |
| Other-molecule core hit | `InterChain` (`BlockKind::InterChain`) |
| Both present | SelfBlocked-wins (`OverlapField::block_kind` returns `SelfBlocked` as soon as one same-mol core hit is seen) |

A public growth run with `CbmcGrow::with_soften_after(2)` must emit at
least one `StepInfo.radscale` `< 1.0`. Every drop between consecutive
rounds is exactly one `0.97` rung, floored at `0.8`.

## Defect (2026-09-04)

Observed 2026-09-04: `radscale = 1.000` for 960 rounds with
`soften_after = 2` because `Chain.deadends` reset on every successful
commit, so `deadends.is_multiple_of(soften_after)` never tripped through
interleaved successes. The public contract after the fix is the
cumulative watermark above, not a consecutive streak.

## Owning tests

Do not copy these. Re-run them to refresh the pin:

- `grow_ladder_fires_after_intermittent_commits` — ladder fires with
  `CbmcGrow::with_soften_after(2)`; `StepInfo.radscale` drops by `0.97`
  rungs and stays ≥ `0.8`.
- `field_block_kind_self_vs_inter` — `OverlapField::block_kind`
  `SelfBlocked` / `InterChain` / `None` / SelfBlocked-wins.

```text
cargo test -p molcrafts-molpack --lib --tests -- grow_ladder_fires_after_intermittent_commits
cargo test -p molcrafts-molpack --lib --tests -- field_block_kind_self_vs_inter
```
