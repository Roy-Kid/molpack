---
title: special-bonds-05-target — Target::with_special_bonds 每目标豁免表
status: done
created: 2026-09-04
depends_on: [special-bonds-03-sink, special-bonds-04-residual]
chain: special-bonds（05 of 6）
grilled: true
---

# special-bonds-05-target — Target::with_special_bonds 每目标豁免表

## Summary

把分子内豁免表收成每目标数据：`Target::with_special_bonds(self, BondDistanceWeights) -> Self` 写入公开字段 `special_bonds`（默认 ≡ 深度 3，`[0,0,0,1]`）。生长只接受二值权重；分数项在 `InternalTree::from_frame` 以 `GrowError::NonBinarySpecialBond` 具名拒绝，入口再包成既有的 `PackError::Grow`。连续与格生长都经 `pub(crate) tree_from_target` 读该目标的表。Rust 与 Python 绑定同步删除 `CbmcGrow::with_exclusion_depth`。04 的 `IntraResidual::from_targets` / assemble 改读 `Target.special_bonds`。Python `Target.with_special_bonds` 与文档转向归 06。

## Domain basis

- 表形 Cassandra `[s12,…,s1N]`，默认 `[0,0,0,1]` ≡ 深度 3。LAMMPS / OPLS / TraPPE / Cassandra 出处见 01。
- 分数权重与 `min_hard_scale` 0.8 相乘无统计意义；元素感知属于 `with_atom_radius`。
- Boon 2017：豁免一类对 = 先验已拥有该类距离。一阶先验拥有 k≤3；没有任何已交付先验拥有 1-6。
- 两个缺陷两个修法：(i) 显式氢 `with_atom_radius(H, ≈0.85)`；(ii) syn-pentane → 二阶先验（grow-axes）。实测 2026-09-04：深度 3 + H=0.85 → 142 轮 / 0.3 s。**不得**把 `[0,0,0,0,0,1]` 写成全原子默认。
- CG：KG c∞ = 1.76 来自 1-3 排除体积。夹具写 `Target::with_special_bonds(from_exclusion_depth(2))`。

## Design

**builder 返回 `Self`。** 与半径族一致。分数表可以出现在公开字段上；具名拒绝只发生在生长编译豁免名单时。

**唯一用户可见错误：** `GrowError::NonBinarySpecialBond { index, weight }`，入口 `PackError::Grow { target, source }`。禁止 `PackError::NonBinarySpecialBond`。`binary_violation` 是 `pub(crate)` 谓词，放在 `src/grow/config.rs`。跑它的地方是 `InternalTree::from_frame`。

**Target：** `pub special_bonds: BondDistanceWeights`。没有 `resolved_special_bonds`。`from_parts` 是默认 3 的唯一生产点。没有 `Target::with_exclusion_depth`。crate root `pub use molrs::BondDistanceWeights`，不进 prelude，不 `type SpecialBonds = …`。rustdoc 第一句：这不是 `ForceField::special_bonds`。

**`InternalTree::from_frame(frame, &weights)`** 在编译 exclusions 之前跑 `binary_violation`。`pub(crate) fn tree_from_target` 只做 template + `&target.special_bonds` + `from_frame`。GrowStage 与两处格点都走它。禁止 min-fold。

**删除引擎旋钮，不留别名。** `GrowConfig.exclusion_depth` / `CbmcGrow::with_exclusion_depth` 删除。Python 解钩（不是窗口）：`python/src/entry.rs` 五处、`molpack.pyi`、`test_grow.py` builder 链。不在 05 给 Python Target 加 `with_special_bonds`。

**04 读者：** assemble / `from_targets` 读 `t.special_bonds`，去掉全局表参数。

**DRAFT：** 删除 `dg-refine.md` 的 `intra_exclusion_depth` / `with_intra_overlap(..., exclusion_depth)` **以及** `dg-refine.acceptance.md` ac-003（改成两个 Target 各一张表）。删除 `grow-axes.md` 的 `WalkGrow::with_exclusion_depth`。`chain-growth-solver.md` 一行：表现在在 `Target.special_bonds`。

**P8 两表测试：** `tests/grow.rs` 用 `GrowStage::from_targets` 两个模板、两张不同的表，断言 exclusions 与残差分类不同。不要把「不折叠」挂在单树 `from_frame` 测试上。

**Reuse decision**

- `with_short_radius` — **pattern**（消费 `with_*`、字段在 Target；不套 `resolved_*`）。
- `InternalTree` 深度入口 — **generalize**。
- `OverlapField` `&[u32]` — **reuse**（不改 `field.rs`）。
- `GrowStage::from_targets` — **reuse**。
- `IntraResidual::from_targets` — **generalize** 读 `Target.special_bonds`。
- min-folds — **new — 不在本段修**。

## Files to create or modify

- `src/target.rs`
- `src/lib.rs`
- `src/grow/mod.rs`
- `src/grow/internal.rs`
- `src/grow/driver.rs`
- `src/grow/config.rs`
- `src/grow/entry.rs`
- `src/grow/lattice/mod.rs`
- `src/grow/lattice/entry.rs`
- `src/pipeline/mod.rs`
- `src/entry/result.rs`
- `tests/target.rs`
- `tests/grow.rs`
- `python/src/entry.rs`
- `python/python/molpack/molpack.pyi`
- `python/tests/test_grow.py`
- `.claude/specs/dg-refine.md`
- `.claude/specs/dg-refine.acceptance.md`
- `.claude/specs/grow-axes.md`
- `.claude/specs/chain-growth-solver.md`（一行）
- `regressions/special-bonds-05-target.md` (new)

## Tasks

- [x] Write failing unit tests for `Target::with_special_bonds` (default `[0,0,0,1]`; depth-2 stores `[0,0,1]`; fractional table is stored and builder returns `Self`)
- [x] Implement `Target.special_bonds`, `with_special_bonds -> Self`, `from_parts` default, crate-root `BondDistanceWeights` re-export
- [x] Write failing unit tests for `InternalTree::from_frame` (`NonBinarySpecialBond` at index 2) and for `GrowStage::from_targets` with two templates / two tables (exclusions and residual classes differ; not min-folded)
- [x] Generalize `InternalTree::from_frame` to run `binary_violation`; add `tree_from_target`; wire GrowStage and both lattice sites; point assemble / `IntraResidual::from_targets` at `Target.special_bonds`
- [x] Delete `GrowConfig::exclusion_depth` / `CbmcGrow::with_exclusion_depth`; migrate tests/grow.rs :1189,:2055,:2283
- [x] Unhook Python `CbmcGrow.with_exclusion_depth` in entry.rs, molpack.pyi, test_grow.py
- [x] Delete leftover exclusion-depth knobs in dg-refine.md + dg-refine.acceptance.md ac-003 rewrite, grow-axes.md; one line in chain-growth-solver.md
- [x] Add regression example `regressions/special-bonds-05-target.md`
- [x] Run full check + test suite including `cargo check --manifest-path python/Cargo.toml`

## Testing strategy

- Happy path：默认表 `[0,0,0,1]`；C12 `exclusions(0) == [0,1,2,3]`。
- Edge：`with_special_bonds` 存 0.5；`from_frame` 拒绝 index 2，Display 含 `1-4` 与 `with_atom_radius`。
- KG 夹具旋钮改为 Target 拼写，断言不变。
- 两表 `GrowStage::from_targets` 钉住 P8。

## Out of scope

- Python `Target.with_special_bonds` 与文档（06）。
- 新的 `PackError` 变体、改 `helpers.rs`。
- 改 `field.rs` 签名。
- CLI / `.inp`。二阶先验。实现 DgRefine / WalkGrow。
- 把 `[0,0,0,0,0,1]` 写成全原子默认。
- [x] Hygiene cleanup (`/mol:simplify` — fmt+clippy clean)
