---
title: special-bonds-02-ladder — 修复生长按链软化阶梯，死端分自阻塞与链间
status: code-complete
created: 2026-09-04
depends_on: []
chain: special-bonds（02 of 6；本段无 molrs 依赖）
grilled: true
---

# special-bonds-02-ladder — 修复生长按链软化阶梯，死端分自阻塞与链间

## Summary

生长驱动里按链软化阶梯今天永远走不到：`Chain.deadends` 在每次成功提交时清零，于是 `deadends.is_multiple_of(soften_after)` 在实测中连续 960 轮仍给出 `radscale = 1.000`（即便 `soften_after = 2`）。本 spec 把这一只计数器拆成三个时钟，并把三个读取器写进契约：连续死端只喂回撤深度与地板上的强制放置；累计死端加已挣档数只喂软化水位；累计永不进入指数回撤。死端分成自阻塞与链间两类具名失败；`commit` 在 `Probe::Blocked` 臂、rollback **之前**分类。既有四条阶梯 pin 一字不改仍绿。本段不是 grow-axes B1：硬核尺度仍是全局标量。

## Domain basis

- **软化下限 0.8 是 Auhl 的 0.8σ 推开地板，本 spec 不改它。** Auhl et al., *J. Chem. Phys.* **119**, 12718 (2003), doi:10.1063/1.1628670。
- **指数回撤喂错时钟会把 D-01 整链撤回请回来。** notes.md D-01：`retract·2^(deadends/4)` 在 40 次死路时达 10240 步。连续时钟（提交即清）才是回撤与「此刻楔死」的合法输入。
- **自阻塞 vs 链间是几何事实。** `OverlapField::probe` 已经按 `mol_of[q] == mol` 区分同分子伙伴。D-01 KG 熔体：1-4 自阻塞 4.4%/trial，分子间 76%。
- Packmol 容差只管分子间（Martínez et al. 2009, doi:10.1002/jcc.21224）。PEO 氢半径归 05/06；本段只在 abort `log::warn!` 里点名 `with_atom_radius` 与 special bonds。

## Design

**三个时钟，三个读取器。** 删除 `Chain.deadends`：

| 字段 | 增加 | 清零 | 唯一读取器 |
|---|---|---|---|
| `deadend_streak` | 每次非-force 失败 | 成功提交置 0 | (1) 回撤 (2) 地板上的 force |
| `deadends_total` | 同上 | 永不因提交/挣档清零 | 软化水位，与 `rungs_earned` 一起 |
| `rungs_earned` | 本链领到一档时 +1 | 一次 run 内不减 | 软化水位 |

```
fn retract_depth(base, streak) -> base.saturating_mul(1 << (streak/4).min(12))
fn rung_due(total, rungs_earned, soften_after) -> total >= (rungs_earned+1)*soften_after  // saturating
fn force_due(hard_scale, min, streak, soften_after) -> hard_scale <= min && streak >= 2*soften_after
```

**禁止交叉接线。** `1 <<` 不得出现 `deadends_total`。`force_due` 不得出现 `deadends_total` / `rungs_earned`。`rung_due` 不得出现 `deadend_streak`。

**时序：** (1) 用尚未计入本轮失败的 streak 算 `force_due`；(2) 成功提交只清 streak；(3) force 则 `force_place`，不递增时钟；(4) 否则两只死端 +1，可能挣档，然后用递增后的 streak 回撤。

**`Probe::Blocked` 保持单元变体；probe 仍 first-hit。** `BlockKind { SelfBlocked, InterChain }` + `block_kind(...) -> Option<BlockKind>` 只在整段失败后调用。SelfBlocked 赢：见到 SelfBlocked 立刻返回，InterChain 记住后扫完。

**`score_atoms` 返回 `Result<(Vec<PlacedAtom>, F), DeadEnd>`。** 拒绝原因在拒绝点产生：`restraints.violated` → `DeadEnd::Restraint`；`Probe::Blocked` 臂里立刻 `block_kind`（propose 不 insert，无需 rollback）。`propose -> Result<Proposal, DeadEnd>`：`trials` 空时返回最后一次失败 trial 的 `Err`。禁止第二遍 walk、禁止 out-parameter。`DeadEnd { Overlap(BlockKind), Restraint }`。

**`commit -> Result<(), BlockKind>`。** 在 `Probe::Blocked` 臂、rollback **之前**分类（部分插入仍在场，同分子兄弟是 SelfBlocked）。禁止 rollback 后再 query（清场后 `block_kind` 为 `None`）。禁止 `expect("commit failed on Probe::Blocked")`。

**这不是 grow-axes B1。** `hard_scale` 仍是 `GrowStage::run` 里的一只标量。`min_hard_scale.fold` / `serial.first()` 不拆。Python `with_soften_after` 文档仍写 consecutive，直到 06。

**报告：** abort `log::warn!` 点名自阻塞比例与 `with_atom_radius` / special-bond。不加 PackResult / StageOutcome / StepInfo 字段。

**Reuse decision**

- `Chain.deadends` — **generalize** 为三时钟。
- `Probe::Blocked` — **generalize**：单元变体保留；`BlockKind` 是查询视图。
- `probe` / `nearest` — **reuse**。
- `min_hard_scale.min()` / `serial.first()` — **reuse as-is**（P8，grow-axes）。

## Files to create or modify

- `src/grow/field.rs`
- `src/grow/moves.rs`
- `src/grow/driver.rs`
- `src/grow/config.rs`
- `src/grow/entry.rs`
- `tests/grow.rs`
- `regressions/special-bonds-02-ladder.md` (new)

## Tasks

- [x] Write failing unit tests for `OverlapField::block_kind` (`tests/grow.rs` → `field_block_kind_self_vs_inter`; existing `Probe::Blocked` asserts stay unedited)
- [x] Generalize `Probe::Blocked` in `src/grow/field.rs` with `BlockKind` + `block_kind` (unit variant and first-hit `probe` stay)
- [x] Write failing unit tests for the three clock readers (`src/grow/driver.rs` `#[cfg(test)]`: `retract_depth_reads_deadend_streak`, `force_due_reads_streak_at_floor`, `rung_due_uses_watermark`, `missed_rung_still_due`) and `tests/grow.rs` → `grow_ladder_fires_after_intermittent_commits`
- [x] Generalize `Chain` in `src/grow/moves.rs`; `score_atoms -> Result<(Vec<PlacedAtom>, F), DeadEnd>`; `propose -> Result<Proposal, DeadEnd>`; `commit -> Result<(), BlockKind>` classifying in the `Probe::Blocked` arm before rollback; SelfBlocked-wins across commit candidates; wire the three readers and abort `log::warn!` in `src/grow/driver.rs`
- [x] Add rustdoc (dimensionless `hard_scale`; three readers; `with_soften_after` is cumulative not consecutive; Python docs stay consecutive until 06)
- [x] Add regression example `regressions/special-bonds-02-ladder.md` (public API only; hard-coded goldens)
- [x] Verify against `grow_cg_kremer_grest_c_inf` first, then `assert_radscale_ladder` / `grow_softening_needs_repeated_dead_ends` / `grow_dense_strict_core_terminates_unconverged` unedited
- [x] Run full check + test suite
- [x] Hygiene cleanup (`/mol:simplify` — stale test comments past-tensed; stencil extract / tests/grow.rs split deferred to `/mol:refactor`)

## Testing strategy

公开生长行为在 `tests/grow.rs`；谓词算术在 `src/grow/driver.rs` `#[cfg(test)]`。

- **Happy path：** `block_kind` 自阻塞 / 链间 / 排除表上 → None / 两者都在 → SelfBlocked。`soften_after = 2` 时 `radscale < 1.0` 至少一次，每档 0.97。
- **Edge cases：** `force_due(0.8, 0.8, 0, 50) == false`（streak 0 不得 force）。`retract_depth` 读 streak 不读 total。`DeadEnd::Restraint` 不加两个桶。
- **既有 pin 一字不改：** 先跑 KG。若 `c_n = 1.76 ± 10%` 因阶梯发火而出带：停下来报告，不改断言。
- **Regression：** `regressions/special-bonds-02-ladder.md`。缺陷注释：2026-09-04 `radscale = 1.000` for 960 rounds with `soften_after = 2`。

## Out of scope

- 不改排除深度、Topology、PackResult、molrs。
- 不把 `hard_scale` 局部化到链（grow-axes B1）。
- 不展开 P8 min-folds。
- Python `api-reference.md` 的 consecutive 文案归 06。
- 不修 PEO 氢半径。
