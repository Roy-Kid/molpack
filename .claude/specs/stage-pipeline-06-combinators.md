---
title: stage-pipeline-06-combinators — Invariant、Layers 与两个管线组合子
status: approved
created: 2026-09-02
chain: stage-pipeline（06 of 7）
---

# stage-pipeline-06-combinators

## Summary

给线性管线加两个组合子：`Repeat`（把一段阶段体重复到 `Until` 满足）与 `Guarded`（用一组 `Invariant`
守住某个阶段的出口，违约时按 `OnViolation` 重跑同一阶段或具名失败）。同时引入 `Layers`（L0–L5
不可修复度位集，定义在它唯一消费者 `Invariant::layer()` 旁边）、`Invariant` trait、`Violation`
与内置的 `RestraintsSatisfied`。`Guarded` **只能**重跑同一阶段或具名失败，永不切换到另一个算法——
用户挑方法，molpack 不替他挑（law P8）。

## Design

**层次**：管线层（`src/pipeline/`）+ 一个新的顶层不变量模块 `src/invariant.rs`。

**`src/invariant.rs`**（预算 ≤ 200 行）：

    pub struct Layers(u8);   // Debug + Clone + Copy + PartialEq + Eq
    // 常量 L0_CONNECTIVITY … L5_LOCAL_GEOMETRY + contains / union / name()
    pub trait Invariant: Send {
        fn name(&self) -> &'static str;
        fn layer(&self) -> Layers;
        fn check(&self, state: &PackState) -> Vec<Violation>;
    }
    pub struct Violation { pub atoms: Vec<usize>, pub what: String }   // Debug + Clone
    pub struct RestraintsSatisfied { tolerance: F }                    // Debug + Clone

- **`Layers` 定义在这里**（评审 🔴 的落实）：04 只把 L0–L5 阶梯写成 `src/stage.rs` 的模块散文，
  因为链内没有任何代码分支于它；本 spec 带来它唯一的读者 `Invariant::layer()`，类型就落在读者
  旁边（law § 3：先有扩展点，再有扩展）。`Guarantees` 至今没有、也不获得 `layers` 字段。
  `src/stage.rs` 的散文改成指向 `crate::invariant::Layers` 的模块引用，不复制常量。
- 扩展形状照 `src/restraint/mod.rs:541-543` 的 "direction-3" 规则：`pub trait` + N 个具体
  `pub struct` 实现，无 `Builtin*` 包装、无标签联合——与 `AtomRestraint` / `Region` / `Handler`
  同一套词汇（law § 1）。
- `RestraintsSatisfied::check` 读**共享 objective** 留在状态上的 `frest`（03 的未缩放裁决原语已经
  在每个阶段末尾把它填好），超过容差即产出一条 `Violation`；不自己算第二套约束度量（law P7 的
  一把尺子）。
- `layer()` 的消费者是 `Guarded` 的失败报告。错误载荷（评审 🟡 的落实）——**层名以渲染后的字符串
  携带**，`error` 因此不获得指向 `invariant` / `stage` 的新边：

      InvariantViolated { stage: &'static str, invariant: &'static str,
                          layer: &'static str, atoms: Vec<usize> }

  `Display` 点名"哪个阶段的哪条不变量在哪一层被破坏、涉及哪些原子"，并给出补救办法（提高
  `max`、放宽容差或换阶段）。

**`src/pipeline/combinators.rs`**（预算 ≤ 250 行）：

    pub enum Until { Passes(usize), Converged }          // Debug + Clone + Copy
    pub enum OnViolation { Fail, Rerun { max: usize } }  // Debug + Clone + Copy

- `Pipeline::with_repeat(body: Vec<Box<dyn StageFactory>>, until: Until)`；
  `Pipeline::with_guarded(stage: impl StageFactory + 'static,
   invariants: Vec<Box<dyn Invariant>>, on_violation: OnViolation)`。
  命名沿用 crate 无例外的 `with_*` 消费型 builder 约定；`Vec<Box<dyn StageFactory>>` 直接可用，
  因为 05 已写下 `impl<T: StageFactory + ?Sized> StageFactory for Box<T>`。
- `Repeat` 与 `Guarded` 本身是 `Stage` 的实现（组合子即阶段），因此管线主体不需要为它们分支——
  衔接检查、handler 括号、裁决一律照旧。`Repeat` 的 `guarantees()` 是体内各阶段 `guarantees()`
  的并；`requires()` 是体内第一个阶段的 `requires()`。`Guarded` 直接透传被守卫阶段的两者。
- **可重入性**：两个组合子都建立在 04 写下的 `Stage::run` 可重入契约上——被重复运行的阶段不得
  消耗自己的配置。04 已把 `GencanStage` 的 optimizer 绑定改成借用形（`ResolvedBinding<'a>`），
  本 spec 用一条端到端测试证明它在 `Repeat` 之下成立：`Repeat { Passes(2) }` 包一个带 `ff`
  optimizer 的 GENCAN 阶段，第二遍仍解析出绑定（`ff`-gated，落在 `tests/optimizer.rs`）。
  **另一半前提由 05 给出**：GENCAN 阶段在 `state.placed() == Placed::All` 时从既有放置继续
  而不是重新 `initial()`，所以 `Repeat` 的第二遍确实是**接着**第一遍跑，而不是两次独立打包；
  本 spec 为这一点单立一条数值断言（见 Tasks / ac-002），否则"重复"与"随机重启"在测试里无法区分。
- `Until::Passes(n)`：体运行 n 遍，`softened` 累加；`Until::Converged`：某遍的
  `StageOutcome.converged` 为真即停（其后不再运行）。两者都受 `should_stop` 打断。
- `OnViolation::Rerun { max }`：最多重跑**同一阶段** max 次；仍违约则返回
  `converged == false` 与最后一次的违约记录（不死循环，不换算法）。`OnViolation::Fail`：立刻返回
  `PackError::InvariantViolated`。
- `StepInfo.stage` 在组合子内保持单调：组合子对外报告为**一个**阶段（它的 `name()` 是
  `"repeat"` / `"guarded"`），内层重复不改变 `index` / `total`，否则 04 的单调契约会破。

**earn-complexity 的明账**（law § 3，照 03 的写法立账）：`Repeat` 与 `Guarded` 在本 spec 内除测试
外没有生产调用方。它们的**第一个真实消费者**是已在 packing 分类 rev 2 中裁定的
`dg-refine`——"连接 ⇄ 精修"配方需要 `Repeat { Until::Converged }` 把连接与精修交替到收敛，
而 ring-closure 的守卫退回需要 `Guarded { RingClosed, Rerun { max } }`。本 spec 同时交付第一个
具体不变量 `RestraintsSatisfied`（读状态上的 `frest`），因此它不是空壳扩展点。
**回退条件**：若 `dg-refine` 与 ring-closure 都被取消，`src/pipeline/combinators.rs` 与
`src/invariant.rs` 应整体回退（两个文件 + 两个 `Pipeline` builder + 一个 `PackError` 变体，
无其他模块依赖它们）；这条移除路径记录在此。

### Reuse decision

- `restraint/mod.rs:541-543` 的 direction-3 规则 — **pattern**：`Invariant` 的扩展形状照抄。
- `error.rs:5-56 PackError` — **pattern**：`InvariantViolated` 具名字段 + 点名补救的 `Display`；载荷限定为 `&'static str` × 3 + `Vec<usize>`，不引入 `error → invariant/stage` 边。
- `context/pack_state.rs::evaluate_unscaled`（03）与状态上的 `frest` — **reuse**：`RestraintsSatisfied` 读它，不另算约束度量。
- `src/stage.rs` 的 L0–L5 模块散文（04）— **generalize**：本 spec 把它编码成 `Layers` 类型，定义在唯一消费者 `Invariant::layer()` 旁边；`stage.rs` 的散文改为指向 `crate::invariant::Layers` 的模块引用。
- `stage.rs::Stage`（04）— **reuse**：组合子是 `Stage` 的实现，不是平行抽象；可重入契约直接消费。
- `optimizer/mod.rs::ResolvedBinding<'a>`（04）— **reuse**：`Repeat` 之下的绑定保有由它保证，本 spec 只加证据。
- `pipeline/mod.rs::Pipeline` 的 `Placed` 推进与 GENCAN 的 `Placed::All` 续跑（05）— **reuse**：`Repeat` 的"第二遍接着跑"完全由它提供，组合子不自己注入种子。
- `pipeline/engine.rs::{StageFactory, impl StageFactory for Box<T>}`（05）/ `entry/mod.rs:169 with_handler` / `target.rs:387 with_restraint` — **pattern / reuse**：`with_repeat` / `with_guarded` 的 builder 形状与装箱体的类型。
- `region.rs` 的 `And` / `Or` / `Not` 组合子 — **pattern only**：命名不撞、语义无关（谓词组合 vs 阶段组合）。

## Files to create or modify

- `src/invariant.rs` (new)
- `src/pipeline/combinators.rs` (new)
- `src/pipeline/mod.rs`
- `src/stage.rs`
- `src/error.rs`
- `src/lib.rs`
- `tests/invariant.rs` (new)
- `tests/pipeline.rs`
- `tests/optimizer.rs`

## Tasks

- [ ] Write failing tests in `tests/invariant.rs` and `tests/pipeline.rs`：`Layers` 的 `contains` / `union` / `name`；`RestraintsSatisfied` 在可满足 / 不可满足约束下的 `check` 结果与 `layer()`；`Repeat { Passes(2) }` 使体阶段运行两次且 `softened` 求和；**第二遍从第一遍的放置出发**——第二遍后的 `fdist` ≤ 第一遍后的 `fdist`，且第二遍后的 `positions()` 与 `GenCanPack::seeded_from(pass1_result).run(..)`（同 seed、同设置）逐位相同；`Until::Converged` 在首次收敛后不再运行；`Guarded` 在通过时不改变结果、`OnViolation::Fail` 返回具名错误、`Rerun { max: 2 }` 最多重跑 2 次后 `converged == false` 而非死循环；组合子内 `StepInfo.stage` 仍单调；确认 RED
- [ ] Write the failing `ff`-gated test in `tests/optimizer.rs`：`Repeat { Passes(2) }` 包住一个带 optimizer 的 GENCAN 阶段，第二遍仍解析出它的绑定
- [ ] Add `src/invariant.rs`：`Layers`（L0–L5 位集 + 常量 + `contains` / `union` / `name`）、`Invariant`、`Violation`、`RestraintsSatisfied`；在 `src/lib.rs` 注册与重导出，并把 `src/stage.rs` 的 L0–L5 散文改为指向 `crate::invariant::Layers` 的模块引用
- [ ] Add `src/pipeline/combinators.rs` with `Repeat` and `Until` 并实现 `Stage`（`requires` / `guarantees` 的合成规则、`softened` 累加、`should_stop` 让出）
- [ ] Add `Guarded` and `OnViolation` to `src/pipeline/combinators.rs` 与 `PackError::InvariantViolated { stage, invariant, layer, atoms }` to `src/error.rs`（载荷为三个 `&'static str` + `Vec<usize>`；只重跑同一阶段或具名失败，永不切换算法）
- [ ] Add `Pipeline::with_repeat` and `Pipeline::with_guarded` to `src/pipeline/mod.rs` 并补模块文档（组合子即阶段、阶段标识不因内层重复而改变、可重入契约与 05 的 `Placed::All` 续跑这两条依赖、earn-complexity 明账与回退条件）
- [ ] Add regression scenario `invariant_regression_restraints_satisfied_golden` to `tests/invariant.rs`（固定夹具下 `frest` 与 `Violation` 计数的硬编码金标，容差 1e-12）
- [ ] Run full check + test suite（含 `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored` 五例）

## Testing strategy

- 归属：`Layers` / `Invariant` / `RestraintsSatisfied` 的契约归 `tests/invariant.rs`，组合子的契约归
  `tests/pipeline.rs` 的新增段，`ff` 绑定保有归 `tests/optimizer.rs`（law § 11）。单测门：
  `cargo test -p molcrafts-molpack --lib --tests -- invariant` 与 `-- pipeline`。
- Happy path：`Layers` 位运算与层名；`RestraintsSatisfied` 在满足约束的状态上返回空违约；
  `Repeat(Passes(2))` 的计数；`Guarded` 在通过时结果与不加守卫时逐位相同。
- **续跑门**：`Repeat { Passes(2) }` 的第二遍是第一遍的续跑而不是随机重启——`fdist` 单调不增，
  且第二遍的坐标与 `GenCanPack::seeded_from(pass1_result)` 逐位相同（`to_bits()`）。这条直接
  消费 05 的 `Placed::All` 续跑；没有它，`Repeat` 与"两次独立打包"在测试里无法区分。
- Edge cases：`Passes(0)`（体不运行，具名拒绝或空操作，由实现在 rustdoc 里钉死）；
  `Rerun { max: 0 }` 等价于 `Fail`；不可满足约束下 `Rerun { max: 2 }` 恰好重跑 2 次后停；
  组合子嵌套（`Guarded` 包 `Repeat`）时 `StepInfo.stage` 仍单调；`should_stop` 在体内生效时
  组合子立即让出；`ff` 下 `Repeat` 第二遍仍有 optimizer 绑定。
- Regression scenario：`invariant_regression_restraints_satisfied_golden`，硬编码金标；对应
  `type: runtime` 验收项。
- 既有逐位门：`tests/pipeline.rs` 的单阶段等价三条与两条衔接等价（05 落地）、
  `tests/grow.rs` 确定性组、`examples_batch` 五例（组合子不进预设路径，理应不受影响——用它证明）。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- 更多内置不变量（键长完好、密度均匀、缠结度）——各归其消费 spec。
- `OnViolation` 的第三种策略（换算法 / 降级）——永久不做（law P8）。
- `Guarantees` 增加 `layers` 字段——本链无分支于它的代码。
- 并行阶段、DAG、每阶段预算。
- Python 侧暴露 `Invariant` / 组合子：v1 保持 Rust-only 扩展点（07 只镜像 `Pipeline` 与阶段对象）。
