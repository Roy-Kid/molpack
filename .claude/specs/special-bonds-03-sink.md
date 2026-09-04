---
title: special-bonds-03-sink — 删除平行 Topology，生长改读 molrs
status: code-complete
created: 2026-09-04
depends_on: [special-bonds-01-molrs]
chain: special-bonds（03 of 6）
grilled: true
---

# special-bonds-03-sink — 删除平行 Topology，生长改读 molrs

## Summary

molpack 不再拥有一份自己的模板键图：删除 `src/topology.rs` 与 `tests/topology.rs`，生长改为读取 `molrs::Topology`。几何读取 `frame_positions` 落在新的 crate-root 叶子 `src/template.rs`（只依赖 std + molrs Frame，不进 `src/frame.rs` 的 context 环）。`GrowError` 具名：`NoAtomsBlock`、`BondOutOfRange { a, b, n }`、`NoBonds`、`Disconnected`、既有 `TemplateTooSmall` / `RingTemplate`。删除 `InternalTree::from_frame`：唯一构造器 `from_frame_with_weights`；03 里默认 3 的唯一家是 `GrowConfig.exclusion_depth`。公开面不把 `molrs::Topology` 再导出成 `molpack::Topology`。

## Domain basis

本段不引入新物理、不改 `fdist`。Packmol 容差只管分子间（Martínez 2009, doi:10.1002/jcc.21224）。默认表 `[0,0,0,1]` ≡ 深度 3 是 01 的 Cassandra 尾项约定。`BondDistanceWeights` 是 molrs **core**，不是 `molrs::ff`。

## Design

**一条键图，一个家。** 删除 molpack CSR `Topology` / `TopologyError`。crate root 不再 `pub mod topology`，也不 `pub use molrs::Topology`。

**几何叶子 `src/template.rs`（new）。** `pub(crate) fn frame_positions(frame) -> Result<Vec<[F; 3]>, FramePositionsError>`，`FramePositionsError::NoAtomsBlock`。缺 atoms 或缺 x/y/z → `NoAtomsBlock`。单位 Å。不与 panicking `frame_to_coords` 合并。**禁止**把此函数放进 `src/frame.rs`（那一模块 import `PackContext`，会把生长读路径拖进 `context ↔ frame` 环）。`grow/mod.rs` 只通过 `crate::template` 取坐标，不 `use crate::frame`。`internal.rs` / `decorate.rs` 禁止 import `crate::frame` 与 `crate::template`——它们只走下面的助手。

**`pub(crate) fn topology_for_growth(frame) -> Result<(molrs::Topology, Vec<[F; 3]>), GrowError>`** 在 `src/grow/mod.rs`。只做两件事：调用 `molrs::Topology::from_frame`（不在 grow 里 `from_edges` / 重写 CSR），再调 `crate::template::frame_positions`，把帧内容映射成具名 `GrowError`。**拒绝顺序的唯一实现者**（含 `RingTemplate`）：

`NoAtomsBlock → BondOutOfRange → NoBonds → TemplateTooSmall → Disconnected → RingTemplate`

`validate_template` 只调用它，不再自己判环。环或不连通的模板不得构造 `InternalTree`。`bfs_order` 里第二份 `Disconnected` 改为 `debug_assert`。

**`GrowError`：** 删除 `Topology(TopologyError)`。不引入 `MolRs(String)`，不嵌 `MolRsError`。`NoBonds` Display 沿用今日原文。`BondOutOfRange { a, b, n }` 从帧读出，不解析 Display。

**`InternalTree`：** 删除 `from_frame`、`from_frame_with_depth`、`DEFAULT_EXCLUDE_BONDS`。唯一入口 `from_frame_with_weights(frame, &BondDistanceWeights)`：内部 `topology_for_growth` + `topo.exclusions(weights)` 机械转成 `Vec<Vec<u32>>`。禁止抄 `topology.rs` 的 BFS。

调用点：`GrowStage` 用 `from_exclusion_depth(config.exclusion_depth)`；格点两处显式 `&from_exclusion_depth(3)`，rustdoc 写明 lattice-v1 全原子。03 默认 3 的唯一家：`GrowConfig`。05 拆除旋钮并把格点字面量换成 `Target.special_bonds`。

DRAFT `dg-refine` 消费 `molrs::Topology` 直接，永不 import `topology_for_growth`（P7）。

**Reuse decision**

- molpack `src/topology.rs` — **generalize**（删）。
- `InternalTree` 深度入口 — **generalize** 为吃表。
- `src/frame.rs` — **new — 不放 frame_positions**（环）。新叶子 `src/template.rs`。
- 读点 — **reuse**。

## Files to create or modify

- `src/template.rs` (new)
- `src/grow/config.rs`
- `src/grow/mod.rs`
- `src/grow/internal.rs`
- `src/grow/lattice/decorate.rs`
- `src/grow/lattice/mod.rs`
- `src/grow/lattice/entry.rs`
- `src/grow/driver.rs`
- `src/lib.rs`
- `src/topology.rs` — 删除
- `tests/topology.rs` — 删除
- `tests/grow.rs`
- `docs/architecture.md`
- `.claude/notes/conventions.md`
- `.claude/specs/dg-refine.md`（一行：refine 读 `molrs::Topology`，不读 `topology_for_growth`）
- `regressions/special-bonds-03-sink.md` (new)

## Tasks

- [x] Write failing unit tests for `frame_positions` (`src/template.rs` `#[cfg(test)]`: zigzag round-trip; missing atoms / missing z → `NoAtomsBlock`)
- [x] Write failing unit tests for `GrowError` and `InternalTree` (`tests/grow.rs`: `NoBonds` before `TemplateTooSmall`/`Disconnected`; isolated atom and two 3-atom chains → `Disconnected`; no atoms / missing z → `NoAtomsBlock`; OOB → `BondOutOfRange`; shuffled-bond tree still decomposes and skip set includes root; ring template is `RingTemplate` and does not construct a tree)
- [x] Implement `src/template.rs` `pub(crate) frame_positions`; `topology_for_growth` in `src/grow/mod.rs` (only grow caller of `crate::template`; owns the full refusal sequence); `from_frame_with_weights`; delete `from_frame` / `from_frame_with_depth` / `DEFAULT_EXCLUDE_BONDS` / `src/topology.rs` / `tests/topology.rs` / crate-root Topology re-exports
- [x] Add rustdoc (Å; refusal order; lattice-v1 AA note on the literal 3; 03 default-3 lives only on GrowConfig)
- [x] Update `docs/architecture.md`, `.claude/notes/conventions.md`, and one line in `.claude/specs/dg-refine.md`
- [x] Add regression example `regressions/special-bonds-03-sink.md`
- [x] Verify default-table bitwise pins (assertion bodies unedited; constructor sites may switch to `from_frame_with_weights`)
- [x] Run full check + test suite
- [x] Hygiene cleanup (`/mol:simplify` — fmt+clippy clean; bfs_order no longer returns Disconnected)

## Testing strategy

`pub(crate)` 的 `frame_positions` 在 `src/template.rs` in-module。生长策略在 `tests/grow.rs`。不把 `molrs::Topology` 当 molpack 类型来测。

- **Happy path：** C12 `from_frame_with_weights(..., from_exclusion_depth(3))` 的 `exclusions(0) == [0,1,2,3]`（`internal_exclusions_depth` 断言正文不改）。
- **Edge cases：** 2 原子无键 → `NoBonds`；两个三原子链 → `Disconnected`；环 → `RingTemplate` 且无树。
- **Bitwise pins：** 见 03 修订稿列出的 `tests/grow.rs` / `tests/pipeline.rs` 名单。

## Out of scope

- 04 残差、05 Target 表、06 Python、02 阶梯。
- `GrowError::MolRs` / crate-root `MolRsError`。
- 自动化 `MOLRS_GIT_REF`（P3）。
- 抄 01 的 exclusions BFS；从 `molrs::ff` 取表。
