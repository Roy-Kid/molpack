---
title: stage-pipeline-04-stage — Solver → Stage：seam 的词汇与三个实现者
status: approved
created: 2026-09-02
chain: stage-pipeline（04 of 7）
---

# stage-pipeline-04-stage

## Summary

把 `src/solver.rs` 改名为 `src/stage.rs`，把 `Solver` 升格为 `Stage`：除 `run` 之外声明前置条件
`requires()` 与出口保证 `guarantees()`；阶段之间交换的状态统一成 `&mut PackState`（03 的
`PackState` / `Placed` 随之升为 `pub`）。三个实现者同时改名为 `GencanStage` / `GrowStage` /
`LatticeStage`，`grow/lattice/entry.rs` 从 `mod.rs` 拆出。`Stage::run` 写下**可重入契约**，
`GencanStage` 因此改成跨调用保留自己的 optimizer 绑定。`StepInfo` 增加阶段标识并标
`#[non_exhaustive]`，`Handler` 获得默认空实现的 `on_stage_start` / `on_stage_end`。行为逐位不变。

## Design

**层次**：packing seam 及其三个实现者。trait 签名变更是编译原子的——`Solver` 与 `Stage` 并存会立刻
造出同一概念的平行抽象（law § 1），因此三个实现者必须在同一步改完，这也是本 spec 不能再切小的原因。

**seam 词汇**（`src/stage.rs`，预算 ≤ 300 行；只 import `context` / `handler` / `target` / `error`）：

    pub trait Stage: Send {
        fn name(&self) -> &'static str;
        fn requires(&self) -> Requires;
        fn guarantees(&self) -> Guarantees;
        fn run(&mut self, state: &mut PackState, targets: &[Target],
               budget: &Budget, handlers: &mut [Box<dyn Handler>]) -> StageOutcome;
    }
    pub struct Requires   { pub placed: Placed }        // #[non_exhaustive]
    pub struct Guarantees { pub placed: Placed }        // #[non_exhaustive]
    pub struct StageOutcome { pub converged: bool, pub softened: usize }  // #[non_exhaustive]
    pub struct Budget { max_loops, precision }          // 不变

- **没有 `validate` 钩子**（评审 🟡 的落实）：本链没有一个实现者会覆写它——GENCAN 保留
  `validate_targets`，生长保留 `validate_grow_cell`（继续报 `PackError::Grow { source: NoBox }`，
  消息与既有测试不变）。一个没人实现、又只有管线会调的前置钩子，是调用方能忘的步骤加上未挣得的
  复杂度（law § 3 / § 8）。seam 只有 `name` / `requires` / `guarantees` / `run` 四件事。
- **没有 `Layers`，`Guarantees` 也没有 `layers` 字段**（评审 🔴 的落实）：本链里没有任何代码
  分支于它——05 的衔接检查只读 `placed`，06 唯一的读者是 `Invariant::layer()`，那是另一个生产者。
  L0–L5 不可修复度阶梯以**模块散文**留在 `src/stage.rs` 的文档里（阶梯定义 + "阶段只对自己声明
  的层负责"），`Layers` 这个**类型**由 06 定义在 `src/invariant.rs`，紧挨它唯一的消费者。
- **`StageOutcome` 只有 `converged` / `softened`**（评审 🔴 的落实）：今天 `SolveOutcome.fdist` /
  `.frest` 零读者（`src/entry/mod.rs:341-342` 从 `sys.fdist` / `sys.frest` 填 `PackResult`），
  把它们保留并经 handler 公开，等于让求解器自报裁决，正是 law P7 的字面禁令。需要数字的
  handler 从唯一的家读：`on_stage_end(&StageInfo, &StageOutcome, &PackContext)`，与
  `on_finish(&PackContext)`（`handler.rs:107`）同形。
- `Placed` 从 `crate::context` 引入（03 定义），不在此处重定义；本 spec 把 `PackState` / `Placed`
  由 `pub(crate)` 升为 `pub`，因为 seam 的签名已经使它们公开。
- `Requires` 只有 `placed` 一个字段：本链没有第二个前置条件需要检查（law § 3）。
- `Stage::run` 不再收 `sys` + `x` 两件东西：阶段用 `state.rigid_split_mut()` 取
  `(&mut PackContext, &mut RigidView)`。刚体自由度因此正式退出共享签名。

**`Stage::run` 的可重入契约**（评审 🔴 的落实，rustdoc 写死）：

> 一个 `Stage` 可能在**演化中的状态**上被运行**多次**（05 的多阶段管线、06 的 `Repeat` /
> `Guarded`）。实现者**不得消耗自己的配置**：第二次 `run` 必须与第一次拥有同样的能力。
> 唯一允许被消耗的是它每次自己新建的临时工作区。

今天 `src/gencan/solver.rs:176` 的 `resolve_bindings(std::mem::take(&mut self.optimizers), …)`
违反这条：第二遍起 `ff` optimizer 绑定被静默丢弃，跑出一个没有名字的降级结果（law § 10）。
修法是**借用而非取走**，不是克隆（`OptimizerBinding.optimizer: Box<dyn Optimizer>` 不可克隆）：

    pub struct ResolvedBinding<'a> {
        pub select: &'a OptimizeSelect,
        pub type_indices: Vec<usize>,
        pub optimizer: &'a mut dyn Optimizer,
    }
    pub(crate) fn resolve_bindings<'a>(bindings: &'a mut [OptimizerBinding],
                                       type_names: &[Option<String>]) -> Vec<ResolvedBinding<'a>>

`GencanStage.optimizers: Vec<OptimizerBinding>` 因此终生留在阶段上，每次 `run` 重新解析类型下标
（廉价、且对不同 `targets` 天然正确）。爆炸半径已核实为三个文件：`src/optimizer/mod.rs`、
`src/gencan/phases.rs`（`run_phase` / `run_iteration` 的 `ff` 参数变为
`&mut [ResolvedBinding<'_>]`）、`src/gencan/solver.rs`；`benches/` 不使用 `ff`，无命中。
`src/gencan/entry.rs:225` 的 `mem::take` 是**入口**一次性交接（`run(self)` 消耗入口），不在契约
范围内，保持不变。

**三个实现者**（纯机械迁移，算法体不动）：

| 今天 | 之后 | 位置 |
|---|---|---|
| `gencan::solver::GencanSolver` | `GencanStage` | `src/gencan/solver.rs` |
| `grow::driver::GrowthSolver` | `GrowStage` | `src/grow/driver.rs` |
| `grow::lattice::LatticeSolver` | `LatticeStage` | `src/grow/lattice/mod.rs` |

- `requires()`：三者都返回 `Placed::None`（GENCAN 自己会 `initial()`；push-off 由
  `RigidView::is_seeded` 派生，02 已落地）。
- `guarantees()`：三者都返回 `Placed::All`。
- 命名规则来自 crate 自身：seam 的实现者带 seam 后缀（今天 `*Solver` 对 `Solver`），入口保持裸名词
  （`GenCanPack` / `CbmcGrow` / `LatticeGrow` 不改）。不留 `Solver` / `SolveOutcome` 别名。
- `src/grow/lattice/mod.rs` 今天同时装着 `LatticeStage`（`:40-302`）与 `LatticeGrow` 入口
  （`:312-401`），另两个算法都是分文件的；本 spec 顺带把入口拆到 `src/grow/lattice/entry.rs`
  （切点 ~304 行，librarian 已核实干净），与 `gencan/entry.rs` / `grow/entry.rs` 对齐。

**handler**（`src/handler.rs`，今天 580 行，**本 spec 后预算 ≤ 700 行**——超出 200–400 常带但仍在
800 上限内，明账记于此，law § 10）：`StepInfo` 增加 `pub stage: StageInfo { index, total, name }` 并标
`#[non_exhaustive]`——它在 `gencan/phases.rs:145`、`grow/driver.rs:323`、`grow/lattice/mod.rs:212`
三处构造，外部构造者是 semver 风险。`StageInfo` 的字段形状照 `PhaseInfo`（`handler.rs:16`）。
`Handler` 增加 `on_stage_start(&StageInfo)` 与 `on_stage_end(&StageInfo, &StageOutcome, &PackContext)`，
**provided 默认空实现**，放在 `handler.rs:114` 那条 "v2 additions — default no-op, backward
compatible" 横幅之下并新增一条标注本 spec 的横幅——这正是 `on_phase_start` / `on_phase_end` 的
既有先例。本 spec 内只有单阶段，`index = 0` / `total = 1`；两个钩子的**调用**由 05 负责。

**生命周期接线**：`src/entry/mod.rs::run` 内把 `PackContext` + `RigidView` 包成 `PackState` 传给阶段，
结束后 `into_parts()` 取回。`PackEngine::solver()` 的返回类型改为 `Box<dyn Stage>`；`solver()` /
`prepare()` 的**删除**是 05 的事，本 spec 只换币种。

### Reuse decision

- `handler.rs:104,134 on_phase_start/on_phase_end` — **pattern**：`on_stage_start/end` 逐条照抄其 provided 默认 + 横幅注释形式；`PhaseInfo` 是 `StageInfo` 的结构模型；`on_finish(&PackContext)`（`:107`）是 `on_stage_end` 末参的模型。
- `restraint/mod.rs:541-543` 的 "direction-3" 扩展规则 — **pattern**：`Stage` 是 `pub trait` + N 个具体 `pub struct` 实现，无 `Builtin*` 包装、无标签联合。
- `gencan/solver.rs` / `entry/mod.rs` — **pattern**：`*Stage` 的构造与错误处理照它们的现有形状。
- `solver.rs:130 Budget` — **reuse**：原样搬迁。
- `solver.rs:154 SolveOutcome` — **generalize**：改名 `StageOutcome`，**去掉** `fdist` / `frest` 两个死字段（裁决只在状态上）。
- `optimizer/mod.rs:85 ResolvedBinding` / `:91 resolve_bindings` — **generalize**：改为借用形，使 `GencanStage` 跨调用保留绑定。
- `PackState`（03）/ `RigidView`（02）/ `Topology`（01）— **reuse**：本 spec 只消费，不扩展（`PackState` / `Placed` 升 `pub`）。
- `grow/lattice/mod.rs:40-302 / :312-401` — **pattern**：按 `gencan/entry.rs` 的分文件形状拆分。
- `Layers` — **new — 不在本 spec**：唯一消费者是 06 的 `Invariant::layer()`，类型定义随它落在 `src/invariant.rs`；本 spec 只留 L0–L5 的模块散文。

## Files to create or modify

- `src/stage.rs` (new — 由 `src/solver.rs` 改名而来)
- `src/solver.rs` (删除)
- `src/grow/lattice/entry.rs` (new)
- `src/grow/lattice/mod.rs`
- `src/gencan/solver.rs`
- `src/gencan/entry.rs`
- `src/gencan/phases.rs`
- `src/grow/driver.rs`
- `src/grow/entry.rs`
- `src/grow/mod.rs`
- `src/optimizer/mod.rs`
- `src/context/pack_state.rs`
- `src/handler.rs`
- `src/entry/mod.rs`
- `src/lib.rs`
- `tests/stage.rs` (new)
- `tests/gencan.rs`
- `tests/grow.rs`
- `tests/optimizer.rs`

## Tasks

- [ ] Write failing tests in `tests/stage.rs` using FAKE stages only：`Box<dyn Stage>` 对象安全、`requires` / `guarantees` 的取值与沿链合成、`Requires` / `Guarantees` / `StageOutcome` 的 `#[non_exhaustive]` 构造路径、`StepInfo.stage` 三字段可读、`on_stage_start` / `on_stage_end` 默认实现不 panic、同一个假阶段连跑两次仍保有自己的配置（可重入契约）；确认 RED
- [ ] Write the failing `ff`-gated re-entrancy test in `tests/optimizer.rs`：同一个 `GencanStage` 连续 `run` 两次，第二次仍解析出它的 optimizer 绑定（今天在 `src/gencan/solver.rs:176` 的 `mem::take` 下会红）
- [ ] Rename `src/solver.rs` to `src/stage.rs` and add the seam vocabulary：`Stage`（`name` / `requires` / `guarantees` / `run`，无 `validate`）、`Requires`、`Guarantees`、`StageOutcome`（只含 `converged` / `softened`）、`Budget`；更新 `src/lib.rs` 的模块声明与重导出，并把 `src/context/pack_state.rs` 的 `PackState` / `Placed` 升为 `pub`
- [ ] Add module docs to `src/stage.rs`：L0–L5 阶梯散文与"阶段只对自己声明的层负责"、`Stage::run` 的可重入契约、"`fdist` / `frest` 的权威是运行后的状态，`StageOutcome` 不携带裁决"
- [ ] Make optimizer bindings borrow-based in `src/optimizer/mod.rs`（`ResolvedBinding<'a>` + `resolve_bindings(&'a mut [OptimizerBinding], …)`）并跟随 `src/gencan/phases.rs` 的 `run_phase` / `run_iteration` 签名，使 `GencanStage` 不再 `mem::take` 自己的绑定
- [ ] Convert the three implementors to `GencanStage` / `GrowStage` / `LatticeStage` on `&mut PackState`（`src/gencan/solver.rs`、`src/grow/driver.rs`、`src/grow/lattice/mod.rs`，以及 `src/gencan/entry.rs` / `src/grow/entry.rs` / `src/grow/mod.rs` 的类型名跟随），并在 `src/entry/mod.rs::run` 内组装 / `into_parts` 拆解 `PackState`
- [ ] Split `src/grow/lattice/entry.rs` out of `src/grow/lattice/mod.rs`（`LatticeGrow` 入口迁出，`mod.rs` 只留 `LatticeStage`）
- [ ] Add `StepInfo.stage: StageInfo` with `#[non_exhaustive]` and the two `Handler` hooks to `src/handler.rs`（`on_stage_end(&StageInfo, &StageOutcome, &PackContext)`，provided 默认空实现，带标注本 spec 的横幅注释，文件 ≤ 700 行），更新三处 `StepInfo` 构造点
- [ ] Move the per-stage assertions to their owners：`tests/gencan.rs` 断言 `GencanStage::requires/guarantees`，`tests/grow.rs` 断言两个生长阶段的同两项；`tests/stage.rs` 不引导任何真实算法
- [ ] Add regression scenario `stage_regression_fake_chain_outcome_golden` to `tests/stage.rs`（两个假阶段依次跑在同一 `PackState` 上，硬编码 `Placed` 迁移序列、`softened` 之和与 `converged` 的字面量金标）
- [ ] Run full check + test suite（含 `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored` 五例）

## Testing strategy

- 归属：**seam 自己的契约**归 `tests/stage.rs`，并且只用**假阶段**验证——对象安全、
  `requires` / `guarantees` 的合成、`#[non_exhaustive]`、默认钩子、可重入（law § 11）。
  真实算法对 seam 的取值归各自的 owner：`tests/gencan.rs` / `tests/grow.rs`；GENCAN 的坐标金标
  继续留在 `tests/gencan.rs`，本 spec 不复制。单测门：
  `cargo test -p molcrafts-molpack --lib --tests -- stage`。
- 夹具：`PackContext::new(ntotat, nmol, ntype)` + `PackState::new`（与 `tests/geometry_cache.rs:28`
  同法），假阶段只改 `Placed` 与计数器，不做几何。
- Happy path：假阶段链的 `Placed` 推进；`StageOutcome` 的两字段；`StageInfo` 在 `total = 1` 时
  `index == 0`。
- Edge cases：`Box<dyn Stage>` 的对象安全；`#[non_exhaustive]` 使外部字面量构造不再可能
  （编译期注释说明）；同一假阶段连跑两次的配置保有；`ff` 下同一 `GencanStage` 连跑两次仍解析出
  绑定（`tests/optimizer.rs`，`ff`-gated）。
- Regression scenario：`stage_regression_fake_chain_outcome_golden`，纯假阶段、硬编码字面量结果，
  不含任何 GENCAN 坐标金标；对应 `type: runtime` 验收项。
- 既有逐位门：`tests/gencan.rs`、`tests/packer.rs`（`seeded_run_contract`、
  `free_chain_push_off_deterministic`）、`tests/grow.rs` 三个确定性测试、`examples_batch` 五例。
  `tests/*.rs` 只改类型名，断言一律不动。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- `Layers` 类型（06，随 `Invariant::layer()` 落在 `src/invariant.rs`）；本 spec 只留阶梯散文。
- `Stage::validate` —— 在有实现者需要之前不加（law § 3）。
- `StageFactory`、`PackEngine::stages()`、`solver()` / `prepare()` 的删除（05）。
- `Pipeline` 与组合子（05 / 06）；本 spec 里 `total` 恒为 1。
- `on_stage_start` / `on_stage_end` 的**调用**（05 负责；本 spec 只给出钩子）。
- 算法本身：`initial` / `movebad` / `pgencan` / 生长驱动的任何数值行为。
