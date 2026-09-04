# .claude/specs/INDEX.md — molpack feature specs

Status legend: **DRAFT** (under review) → **APPROVED** (ready to implement) → **IN PROGRESS** → **CODE-COMPLETE** (evaluators owed) → **DONE** (closed and deleted).

Add via `/mol:spec <feature description>`. Implement via `/mol:impl <slug>` (chains: `/mol:impl-all <prefix>`). Close via `/mol:close <slug>`, which deletes the spec, its acceptance file, and this entry. Specs are active artifacts — finished ones do not stay here.

## special-bonds (chain — 01 closed in molrs; 05–06 here)

来源：PEO 全原子生长卡住（2026-09-04）。根因是显式氢半径 1.0 Å，不是要把默认豁免加深到 1-6。表形 Cassandra `[s12,…,s1N]`，默认 ≡ 深度 3；生长只接受二值；`PackResult.intra` 在表可配置之前报告分子内残差。全原子一等建议是 `Target::with_atom_radius(H, ≈0.85)`。

- [special-bonds-05-target](special-bonds-05-target.md) — Target::with_special_bonds；生长编译点具名拒绝分数权重；删除引擎 exclusion_depth [approved]
- [special-bonds-06-mirror](special-bonds-06-mirror.md) — Python Target.with_special_bonds + IntraResidual 镜像；文档把氢半径写成全原子一等建议 [approved]

01 在 molrs 已关闭（`BondDistanceWeights` + `Topology::from_frame` / `exclusions`）。02-ladder 已关闭（三时钟 + `BlockKind`）。03-sink 已关闭（molpack Topology 删除；`topology_for_growth` + `template.rs`）。04-residual 已关闭（PackResult.intra）。04-residual 已关闭（PackResult.intra）。

## packing-taxonomy (chain — `stage-pipeline` 01–07 landed 2026-09-03; next: `grow-axes` / `dg-refine`)

来源：grow 族算法评审 + packing 分类 rev 2（2026-09-02，用户裁定：所有 packer 共享 Stage/Pipeline trait 管理生命周期与 handler；全部化学、力场无关；键长不必精确，后接力场 minimize；生成族按六条正交轴组合）。顺序（architect 裁定）：`stage-pipeline`（含共享叶子 `src/topology.rs`）先落地；`grow-axes` 可在其前后独立落地（预设写在今天的 `PackEngine` 上）；`dg-refine` 前置 `src/objective.rs` 的行为保持拆分，该拆分与 DRAFT `pair-loop-context-split` 合并为一次 objective 重组。

- stage-pipeline（父 spec 已按 large-spec-split 拆为 7 段链并删除，2026-09-02；设计理由分布在各子 spec 的 Design 里；三轮 architect design-mode 记录见本会话）：
  - stage-pipeline-01-topology — DONE 2026-09-03（6/6 verified；`src/topology.rs` 叶子落地，grow 三处消费者改接，`GrowError::Topology` 收敛；随 squash 提交关闭并删除）
  - stage-pipeline-02-view — DONE 2026-09-03（7/7 verified；`src/context/rigid_view.rs` `RigidView { x, nmol }` 吸收 `PlacementsMut` / `init_xcart_from_x` 两个方向 / 种子注入 / 两处生长写回；连续驱动 abort 路径先同步 `xcart`，`grow_abort_writeback_golden` 逐位金标前后皆绿；`Solver::solve` 收 `&mut RigidView`；`push_off` 未动；`initial()` 的临时视图过桥记为 D-05 归 05；随本段提交关闭并删除）
  - stage-pipeline-03-state — DONE 2026-09-03（8/8 verified；`src/context/pack_state.rs` `PackState { ctx, placed, rigid }` / `Placed` 均 `pub(crate)`（测试为 crate 内 `pack_state/tests.rs`，`Debug` 手写）；未缩放裁决原语合一：`evaluate_unscaled` 唯一实现在 context 层，对称存还 `scale`/`scale2` 与 `radius`，三处调用点改接，`molpack::gencan::phases::evaluate_unscaled` 公开路径撤下；随本段提交关闭并删除）
  - stage-pipeline-04-stage — DONE 2026-09-03（10/10 verified；`src/solver.rs` → `src/stage.rs`，`Stage { name, requires, guarantees, run }` 收 `&mut PackState`，`StageOutcome { converged, softened }` 不带裁决；`GencanStage` / `GrowStage` / `LatticeStage`，`grow/lattice/entry.rs` 拆出；optimizer 绑定改借用形消除 `mem::take`（`ff` 再入测试）；`StepInfo.stage` + `on_stage_start/end` 钩子；`PackState`/`Placed` 升 `pub`；docs/ 跟随改名；随本段提交关闭并删除）
  - stage-pipeline-05-pipeline — DONE 2026-09-03（12/12 verified，ac-010 的 tox 门受 D-04 限制留待 CI；`src/pipeline/{mod,engine,bracket}.rs`：`StageFactory` / `PackEngine`（`run` 必需）/ `EngineSetup`（含 `settings` 引用）迁入，`Pipeline::{new, single, with_stage}`，生命周期体 validate / resolve_stages / run_stages / assemble，衔接检查与设置检查在任何 handler 通知前，handler 采纳，`StageTagger` 改写 `StepInfo.stage`，push-off 状态化（`GencanSettings.push_off` 删除），`prepare()`/`solver()` 删除，三处网格前奏共用 `install_resolved_cell`（radmax 从 `radius_ini`）；两组逐位等价门绿；docs 三页跟随；随本段提交关闭并删除）
  - stage-pipeline-06-combinators — DONE 2026-09-03（9/9 verified；`Stage::run` 改为可失败 `Result<StageOutcome, PackError>`（对 04 的修正）；`src/invariant.rs`：`Layers` L0–L5 位集落在唯一消费者 `Invariant::layer()` 旁、`Violation`、`RestraintsSatisfied`（读共享 `frest`）；`src/pipeline/combinators.rs`：`Repeat`/`Until`、`Guarded`/`OnViolation`（同阶段重跑或具名失败，永不换算法），`with_repeat`/`with_guarded` 经 `with_stage` 同一采纳路径；`PackError::InvariantViolated`；`Repeat` 第二遍逐位 ≡ `seeded_from`，`fdist` 单调断言撤回；预算修正 combinators ≤ 320 / mod ≤ 500；随本段提交关闭并删除）
  - stage-pipeline-07-bindings — DONE 2026-09-03（7/7 verified；Python `Pipeline([stage, …])` / `.with_stage` / 共享 `with_*` / `run`，每入口 `IntoStageFactory`（`to_stage_factory`）+ 唯一 `stage_entry_registry!`（派发 + `TypeError` 文案），阶段对象的 handler 经 `take_handlers` 被采纳，`StepInfo.stage` → `StageInfo{index,total,name}`；`.pyi`/`_protocols`/`__init__` 同步；五个文档页 + CLAUDE.md + conventions + `/mol:map` 蓝图刷新；Python 门 `tox -c python -e py` 167 过（含 9 条 `test_pipeline.py`，两金标）；顺带修 tox `commands_pre[3]` 路径 bug；随本段提交关闭并删除）
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
