---
title: special-bonds-04-residual — PackResult 报告分子内残差
status: done
created: 2026-09-04
depends_on: [special-bonds-03-sink]
chain: special-bonds（04 of 6）
grilled: true
---

# special-bonds-04-residual — PackResult 报告分子内残差

## Summary

每一次引擎运行的 `PackResult` 在组装时附上分子内残差：同链**计分对**的最小距离，以及同链**豁免对**的最小距离（Å，最小映像）。分类尺是调用方传入的 `BondDistanceWeights`，不是钉在构造函数里的常数；`Pipeline::assemble` 在本段对每个 target 传入 `from_exclusion_depth(3)`。05 只改 assemble 的列表构造为 `Target.special_bonds`，不改 `IntraResidual` 的含义。`fdist` 仍然只看分子间。Python 镜像归 06。

## Domain basis

- Packmol 容差只管分子间（Martínez 2009, doi:10.1002/jcc.21224）。`src/objective.rs:245-246` 同分子 `return None`，所以 `fdist` 对分子内接触恒为 0。
- 豁免不是统计中性的。孤立 PEO C_n：无核 4.77 / 深度 5 为 6.14 / 深度 3–4 幸存 0.1% 为 11.30。Auhl 2003 (doi:10.1063/1.1628670) 要把原子级链的分子内排除体积延伸得更远。
- 几何：1-6 可到 1.71 Å（Amber PEO 模板）。`fdist` 为 0 而分子内 1-6 为 1.71，正是本报告要关掉的静默债。
- 单位：残差 Å；键距是无量纲图距离。

## Design

```
pub struct IntraResidual { pub scored: F, pub exempted: F }  // +∞ if empty; no Default
impl IntraResidual {
    pub fn from_targets(
        targets: &[Target],
        positions: &[[F; 3]],
        simbox: &SimBox,
        tables: &[BondDistanceWeights],
    ) -> Self
}
PackResult { ..., pub intra: IntraResidual }
```

- `from_targets` 不默认表。`debug_assert` 不够：`assert_eq!(targets.len(), tables.len())`（release 也查）。禁止无表重载。
- assemble 调用点：`vec![BondDistanceWeights::from_exclusion_depth(3); n]`。`PackResult.intra` 的 rustdoc 必须写明：04 广播深度 3；05 把列表构造换成 `Target.special_bonds`。不读 `GrowConfig::exclusion_depth`。一段深度 2 的 `CbmcGrow` 与深度 3 报告之间的窗口由 05 关掉。
- rustdoc 对 03/05 的类比用 `InternalTree::from_frame_with_weights` 和 `Target.special_bonds`，**不要**写 `resolved_special_bonds` 或 `from_frame_with_depth`。
- 坐标：target 声明序。盒子：`sys.simbox`。不可失败。
- 身份豁免（所有 i≠j 计分）只用于缺 template 与零边图。`NotFound` / `Validation` 省略该 target（两个 +∞），rustdoc 具名；不把坏掉的 1-2 记成 scored。
- 算法：每副本 `i < j`，`SimBox::shortest_vector_impl`。`j ∈ exclusions[i]` → exempted，否则 scored。**尾项规则：** 不在 skip list 里时，若 `weight(usize::MAX) != 0` 则 scored，否则 exempted（零尾项下远对仍豁免）。禁止 NeighborList / AABB / PackContext grid / OverlapField。
- python/src/result.rs 一字不改。crate root 再导出 `IntraResidual`，不进 prelude。

**Reuse decision**

- `PackResult` — **reuse**（加字段）。
- objective 同分子跳过 — **reuse**。
- `Topology::exclusions` — **reuse**。
- `from_exclusion_depth` — **reuse** 在 assemble 调用点。
- `SimBox::shortest_vector_impl` — **reuse**。
- NeighborList / OverlapField / PackContext grid / ViolationMetrics 字段 — **new — 不用**。

## Files to create or modify

- `src/entry/result.rs`
- `src/entry/mod.rs`
- `src/pipeline/mod.rs`
- `src/lib.rs`
- `src/objective.rs`
- `regressions/special-bonds-04-residual.md` (new)

## Tasks

- [x] Write failing unit tests for `IntraResidual` (`src/entry/result.rs` `#[cfg(test)]`: linear hexamer 1.0/4.0; folded 3-4-5 scored 3.0; five-bead depth-1 vs 3 scored 2.0 vs 4.0; bonded dimer scored +∞; coincident unbonded scored 0; second copy ignored; PBC scored 0.5; bondless dimer scores the pair; Validation omit not identity)
- [x] Implement `IntraResidual::from_targets` (`tables: &[BondDistanceWeights]`, skip-list + tail rule, one MIC `i<j` loop)
- [x] Add `PackResult.intra` filled in `Pipeline::assemble` from `vec![from_exclusion_depth(3); n]` + `&sys.simbox` + target-order positions
- [x] Re-export `IntraResidual` from `src/entry/mod.rs` and `src/lib.rs` (not prelude)
- [x] Add rustdoc (Å; table argument; assemble broadcasts depth 3; 05 swaps the list builder for `Target.special_bonds`; named Topology errors omit)
- [x] Add regression example `regressions/special-bonds-04-residual.md`
- [x] Verify linear-hexamer, folded 3-4-5, and five-bead depth-1 vs 3 goldens
- [x] Run full check + test suite
- [x] Hygiene cleanup (`/mol:simplify` — fmt+clippy clean)

## Testing strategy

in-module 测试只打 `from_targets`。手算金标见修订稿：线性六聚体 1.0/4.0；折叠 3-4-5 的 1-6 = 3.0；五珠深度 1 vs 3 钉住「构造函数不是常数」。

## Out of scope

- 不读 `Target.special_bonds`（05）。不读引擎 `exclusion_depth`。
- 不改 Python、阶梯、`fdist` 跳过、半径。
- 不引入 NeighborList。
