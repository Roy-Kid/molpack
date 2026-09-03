# .claude/specs/INDEX.md — molpack feature specs

Status legend: **DRAFT** (under review) → **APPROVED** (ready to implement) → **IN PROGRESS** → **CODE-COMPLETE** (evaluators owed) → **DONE** (closed and deleted).

Add via `/mol:spec <feature description>`. Implement via `/mol:impl <slug>` (chains: `/mol:impl-all <prefix>`). Close via `/mol:close <slug>`, which deletes the spec, its acceptance file, and this entry. Specs are active artifacts — finished ones do not stay here.

## packing-taxonomy (chain — implement in order, start with `stage-pipeline`)

来源：grow 族算法评审 + packing 分类 rev 2（2026-09-02，用户裁定：所有 packer 共享 Stage/Pipeline trait 管理生命周期与 handler；全部化学、力场无关；键长不必精确，后接力场 minimize；生成族按六条正交轴组合）。顺序（architect 裁定）：`stage-pipeline`（含共享叶子 `src/topology.rs`）先落地；`grow-axes` 可在其前后独立落地（预设写在今天的 `PackEngine` 上）；`dg-refine` 前置 `src/objective.rs` 的行为保持拆分，该拆分与 DRAFT `pair-loop-context-split` 合并为一次 objective 重组。

- stage-pipeline（父 spec 已按 large-spec-split 拆为 7 段链并删除，2026-09-02；设计理由分布在各子 spec 的 Design 里；三轮 architect design-mode 记录见本会话）：
  - stage-pipeline-01-topology — DONE 2026-09-03（6/6 verified；`src/topology.rs` 叶子落地，grow 三处消费者改接，`GrowError::Topology` 收敛；随 squash 提交关闭并删除）
  - [stage-pipeline-02-view](./stage-pipeline-02-view.md) — `src/context/rigid_view.rs` `RigidView { x, nmol }` 吸收 `PlacementsMut`、`init_xcart_from_x` 两个方向、种子注入、两处生长写回；连续驱动 abort 路径先把 `xcart` 变成唯一的家；`push_off` 不动 — **APPROVED**
  - [stage-pipeline-03-state](./stage-pipeline-03-state.md) — `src/context/pack_state.rs` `PackState { ctx, placed, rigid }`（`pub(crate)`，包裹不抽取）；未缩放裁决原语的唯一归属（三处副本合一，`scale`/`scale2` 对称存还） — **APPROVED**
  - [stage-pipeline-04-stage](./stage-pipeline-04-stage.md) — `src/solver.rs` → `src/stage.rs`：`Stage { name, requires, guarantees, run }`、`StageOutcome { converged, softened }`；三个 `*Stage` 实现者；`Stage::run` 可重入契约（`GencanStage` 借用形 optimizer 绑定）；`StepInfo.stage` + 两个 handler 钩子 — **APPROVED**
  - [stage-pipeline-05-pipeline](./stage-pipeline-05-pipeline.md) — `src/pipeline/{engine,mod}.rs`：`StageFactory`（含 `take_handlers`）/ `PackEngine`（`run` 必需）/ `EngineSetup` 迁入，`Pipeline::{new, single, with_stage}`，衔接检查、handler 采纳、settings 具名拒绝、push-off 状态化（`Placed::All`）、阶段边界缓存失效；两组 `≡ seeded_from` 逐位门 — **APPROVED**
  - [stage-pipeline-06-combinators](./stage-pipeline-06-combinators.md) — `src/invariant.rs`（`Layers`、`Invariant`、`RestraintsSatisfied`）+ `src/pipeline/combinators.rs`（`Repeat`/`Until`、`Guarded`/`OnViolation`，永不换算法） — **APPROVED**
  - [stage-pipeline-07-bindings](./stage-pipeline-07-bindings.md) — Python `Pipeline` 镜像（每入口 `IntoStageFactory` + 单一注册表）、`StepInfo.stage`、`.pyi`/`_protocols`、五个文档页 + 仓库知识 + `/mol:map` — **APPROVED**
- [dg-refine](./dg-refine.md) — 笛卡尔距离几何精修阶段 `DgRefine`：软核 overlap 项 + 分子内重叠项 + 1-2/1-3 键距弹簧 + 可选手性/约束项，复用 gencan 线搜索原语与 `gxcar`，半径阶梯 s₀→1；熔体密度 push-off 不再出境到 MD — **DRAFT**
- [grow-axes](./grow-axes.md) — 一个 `Grow<Space, ExcludedVolume>` 驱动 + Selector / Escape / Schedule 三个 enum；CbmcGrow / LatticeGrow 成为预设（后者获得哈希流与轮转），新增 `WalkGrow`（理想链，Auhl 路线第一步）；删除 relax / serial / void_bias 旗标，软化按链局部化 — **DRAFT**

## Active

- [pair-loop-context-split](./pair-loop-context-split.md) — move the pair kernel's read set into a `PairInputs` sub-struct so the accumulators stay borrowable, then share the five duplicated traversals. Architecture first; perf recorded per step on a dedicated node, not gated. Now paired with the `objective.rs` split that `dg-refine` needs. **DRAFT**
- [chain-growth-solver](./chain-growth-solver.md) — a `Solver` seam in `pack_with_report`, and a configurational-bias chain-growth solver ranked alongside `gencan` on it: per-target method selection (`Target::with_method`), the box at final volume from step 0, chains grown a torsion at a time under a mandatory geometric torsion prior (RIS/C∞-calibrated — uniform sampling is quantitatively wrong per litrev), recoil-ready retraction, AA + CG. No `ff` dependency. **APPROVED** (revised 2026-08-28: litrev + four design principles; Tasks 1–11, 13 landed — Task 12 measurement and final acceptance ledger outstanding)
- [lattice-growth-phase](./lattice-growth-phase.md) — 金刚石格相生长:格上 SAW 构造 + RIS-MC 修复 + 装饰回连续,melt+ 密度秒级 — v1 LANDED(线性全管线);分支/MC修复/统计验收仍 DRAFT（MC 修复现归精修族 R3，见 grow-axes / packing-taxonomy）
- [collective-com-restraints](./collective-com-restraints.md) — packing to collective targets (COM-level distribution restraints) — **DRAFT**
- [triclinic-cell-downshift](./triclinic-cell-downshift.md) — pack into arbitrary lattices (true triclinic minimum image; lifts the growth entries' orthorhombic-only restriction) — **DRAFT**

## Landed — pending `/mol:close`

Files still on disk; `/mol:close <slug>` verifies the ledger and deletes them.

- [engine-entry-split](./engine-entry-split.md) — 删除大一统 Molpack:按算法独立入口(GenCanPack/CbmcGrow/LatticeGrow)共享生命周期 trait;拓扑一等;无历史包袱 — DONE 2026-09-01(含 placement-seeding 收尾;examples_batch release 门通过)。Its "no Pipeline" Non-goal is superseded by `stage-pipeline`.
- [placement-seeding](./placement-seeding.md) — PackResult 携带放置解,GenCanPack::seeded_from 跨入口逐位接续;门槛 2 原裁定落地,with_push_off/after_solve 退役 — DONE
- [gencan-anneal-budget](./gencan-anneal-budget.md) — pack_solvprotein 5.5× slowdown: real root cause was a missing `avoid_overlap` fixed-atom rejection in `initial.rs` (solvent seeded inside the fixed protein); fixed → now faster than packmol. Proposed anneal-budget workaround SUPERSEDED — RESOLVED
