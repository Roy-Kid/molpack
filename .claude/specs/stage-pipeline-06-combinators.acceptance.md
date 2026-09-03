---
slug: stage-pipeline-06-combinators
criteria:
  - id: ac-001
    summary: Layers lands next to its only consumer; the combinators have their own file
    type: code
    pass_when: |
      `src/invariant.rs` 定义 `pub struct Layers`（`Debug + Clone + Copy + PartialEq + Eq`，
      L0–L5 常量 + `contains` / `union` / `name`）、`pub trait Invariant`（`name` / `layer` /
      `check`）、`pub struct Violation`、`pub struct RestraintsSatisfied`，文件 ≤ 200 行；
      `grep -n 'struct Layers\|Layers(' src/stage.rs` 无命中——`src/stage.rs` 只允许出现指向
      `crate::invariant::Layers` 的文档引用，不得定义或构造该类型；
      `src/pipeline/combinators.rs` 定义 `Repeat` / `Until` / `Guarded` / `OnViolation`
      且 `Repeat` 与 `Guarded` 均实现 `Stage`，文件 ≤ 250 行；
      `grep -rn 'use crate::\(grow\|gencan\|refine\|initial\)' src/invariant.rs src/pipeline/` 无命中。
    status: pending

  - id: ac-002
    summary: Repeat runs the body, accumulates honestly, and continues from the previous pass
    type: runtime
    pass_when: |
      `Repeat { until: Until::Passes(2) }` 使体阶段运行两次（handler 计数为 2）且
      `PackResult.softened` 等于两次之和；`Until::Converged` 在首次 `converged` 之后不再
      运行体（handler 计数为 1）；`Passes(0)` 与 `Rerun { max: 0 }` 的行为与 rustdoc 一致；
      且第二遍是**续跑**：包住一个 GENCAN 阶段时，第二遍后的 `fdist` ≤ 第一遍后的 `fdist`，
      并且第二遍后的 `positions()` 与 `GenCanPack::seeded_from(&pass1_result).run(..)`
      （同 seed、同共享设置、声明了周期盒的夹具，与 05 ac-008 同一比较口径）逐位相同（`to_bits()`）。
    status: pending

  - id: ac-003
    summary: Guarded reruns the same stage or fails by name — never switches algorithm
    type: runtime
    pass_when: |
      在可满足约束下 `Guarded(stage, [RestraintsSatisfied], _)` 的结果与不加守卫时逐位相同；
      在不可满足约束下 `OnViolation::Fail` 返回 `PackError::InvariantViolated { stage,
      invariant, layer, atoms }`（三个字符串字段 + `Vec<usize>`，`Display` 点名四者与补救
      办法），`OnViolation::Rerun { max: 2 }` 恰好重跑同一阶段 2 次后返回
      `converged == false` 而非死循环或换算法。
    status: pending

  - id: ac-004
    summary: error.rs gains no edge to invariant or stage
    type: code
    pass_when: |
      `grep -n 'use crate::' src/error.rs` 不含 `invariant` / `stage` / `pipeline` / `context`；
      `InvariantViolated` 的 `layer` 字段类型是 `&'static str`（渲染后的层名），不是 `Layers`。
    status: pending

  - id: ac-005
    summary: A repeated stage keeps its ff bindings on the second pass
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests --features ff -- optimizer` 中的新增测试
      通过：`Repeat { Passes(2) }` 包住一个带 optimizer 的 GENCAN 阶段，第二遍仍解析出它的
      绑定（依赖 04 的可重入契约与借用形 `ResolvedBinding`）。
    status: pending

  - id: ac-006
    summary: Stage identity stays monotone through combinators
    type: runtime
    pass_when: |
      在 `Guarded(Repeat(body))` 的嵌套管线下，`StepInfo.stage.index` 仍单调递增且
      `total` 等于管线的顶层阶段数（内层重复不改变两者）。
    status: pending

  - id: ac-007
    summary: The complexity ledger names the first consumer and the rollback
    type: docs
    pass_when: |
      `src/pipeline/mod.rs` 的模块文档写明组合子的第一个真实消费者（`dg-refine` 的
      "连接 ⇄ 精修"配方与 ring-closure 的守卫退回）与回退条件（两个文件 + 两个 builder +
      一个 `PackError` 变体整体移除，无其他模块依赖），并点名组合子依赖的两条前提
      （04 的可重入契约、05 的 `Placed::All` 续跑）；`cargo doc -p molcrafts-molpack
      --no-deps` 零警告。
    status: pending

  - id: ac-008
    summary: The suite and the Packmol regression stay green
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` 通过（既有 RED
      `grow_cg_kremer_grest_c_inf` 除外，未被修改或跳过）；
      `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored`
      五例通过；`tests/pipeline.rs` 的五条等价断言（05 落地）未改且全绿。
    status: pending

  - id: ac-009
    summary: Regression scenario reproduces the hard-coded invariant goldens
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- invariant_regression_restraints_satisfied_golden`
      通过：固定夹具下 `RestraintsSatisfied::check` 返回的违约条数与状态 `frest` 与测试内
      硬编码字面量在 1e-12 内相等；无第三方运行时。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-004** 是结构门：`Layers` 落在唯一消费者旁（`stage.rs` 只留文档引用）、文件预算、
  无向下依赖算法模块、错误层不获得新边。
- **ac-002 / ac-003** 是行为门，其中 ac-002 的续跑断言把"重复"与"随机重启"区分开，ac-003 直接
  编码 law P8（不静默换算法）。
- **ac-005** 证明 04 的可重入契约在组合子之下真的成立。
- **ac-006** 保护 04 建立的阶段标识单调契约。
- **ac-007** 是 earn-complexity 的明账（law § 3 / § 10）。
- **ac-008 / ac-009** 是回归与本 spec 的回归场景。
