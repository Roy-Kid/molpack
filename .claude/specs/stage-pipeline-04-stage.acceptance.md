---
slug: stage-pipeline-04-stage
criteria:
  - id: ac-001
    summary: The seam is Stage with four methods, and Solver is gone without aliases
    type: code
    pass_when: |
      `src/stage.rs` 定义 `pub trait Stage`，方法恰为 `name` / `requires` / `guarantees` /
      `run`（`grep -n 'fn validate' src/stage.rs` 无命中），以及 `Requires` / `Guarantees` /
      `StageOutcome` / `Budget`，文件 ≤ 300 行；`src/solver.rs` 不存在；
      `grep -rn '\bSolver\b\|SolveOutcome' src/ python/src/ docs/` 无命中；
      `grep -n 'Layers' src/stage.rs` 无命中（该类型属于 06）。
    status: pending

  - id: ac-002
    summary: stage.rs imports no algorithm module
    type: code
    pass_when: |
      `grep -n 'use crate::' src/stage.rs` 只出现 `context` / `handler` / `target` /
      `error`；`grep -rn 'use crate::\(grow\|gencan\|refine\|initial\)' src/stage.rs` 无命中。
    status: pending

  - id: ac-003
    summary: StageOutcome carries no verdict; handlers read the state
    type: code
    pass_when: |
      `StageOutcome` 的字段恰为 `converged: bool` 与 `softened: usize`
      （`grep -n 'fdist\|frest' src/stage.rs` 无命中）；`src/handler.rs` 的
      `on_stage_end` 签名为 `(&mut self, &StageInfo, &StageOutcome, &PackContext)`；
      `src/entry/mod.rs` 仍从 `sys.fdist` / `sys.frest` 填 `PackResult`。
    status: pending

  - id: ac-004
    summary: Three implementors renamed; the lattice entry has its own file
    type: code
    pass_when: |
      `GencanStage` / `GrowStage` / `LatticeStage` 分别定义于 `src/gencan/solver.rs` /
      `src/grow/driver.rs` / `src/grow/lattice/mod.rs`，三者 `requires()` 均为
      `Placed::None`、`guarantees()` 均为 `Placed::All`；`src/grow/lattice/entry.rs` 存在并
      持有 `LatticeGrow`；`src/grow/lattice/mod.rs` 不再定义 `LatticeGrow`；
      `PackState` / `Placed` 在 `src/context/pack_state.rs` 已升为 `pub` 并从 `src/lib.rs` 重导出。
    status: pending

  - id: ac-005
    summary: A stage keeps its configuration across runs
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests --features ff -- optimizer` 中新增的
      再入测试通过：同一个 `GencanStage` 连续 `run` 两次，第二次仍解析出它的 optimizer 绑定
      （该测试在本 spec 之前的 `src/gencan/solver.rs:176` `mem::take` 下会红，RED 记录在提交
      信息里）；`grep -n 'mem::take(&mut self.optimizers)' src/gencan/solver.rs` 无命中。
    status: pending

  - id: ac-006
    summary: StepInfo carries the stage and is non_exhaustive; handler.rs stays within budget
    type: code
    pass_when: |
      `src/handler.rs` 的 `StepInfo` 标 `#[non_exhaustive]` 且含
      `pub stage: StageInfo { index, total, name }`；两个 provided 默认空实现的钩子位于标注
      本 spec 的横幅注释之下；三处 `StepInfo` 构造点均已填 `stage`；
      `src/handler.rs` ≤ 700 行。
    status: pending

  - id: ac-007
    summary: The module docs pin the ladder, the re-entrancy contract and the verdict authority
    type: docs
    pass_when: |
      `src/stage.rs` 的模块文档写明 L0–L5 阶梯与"阶段只对自己声明的层负责"（散文，不是类型）、
      `Stage::run` 的可重入契约（"可能被多次运行，不得消耗自己的配置"）、以及
      "`fdist` / `frest` 的权威是运行后的状态"；`cargo doc -p molcrafts-molpack --no-deps` 零警告。
    status: pending

  - id: ac-008
    summary: The seam test file boots no real algorithm
    type: code
    pass_when: |
      `grep -n 'GenCanPack\|CbmcGrow\|LatticeGrow\|GencanStage\|GrowStage\|LatticeStage' tests/stage.rs`
      无命中（`tests/stage.rs` 只用假阶段）；`tests/gencan.rs` 与 `tests/grow.rs` 各自持有
      对应阶段的 `requires()` / `guarantees()` 断言。
    status: pending

  - id: ac-009
    summary: Behaviour is bitwise unchanged across the rename
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` 通过（既有 RED
      `grow_cg_kremer_grest_c_inf` 除外，未被修改或跳过），`tests/*.rs` 中除类型名外
      无断言改动；`cargo test -p molcrafts-molpack --release --features io --test
      examples_batch -- --ignored` 五例通过；`seeded_run_contract` 与
      `free_chain_push_off_deterministic` 全绿。
    status: pending

  - id: ac-010
    summary: Regression scenario reproduces the hard-coded fake-stage goldens
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- stage_regression_fake_chain_outcome_golden`
      通过：两个假阶段依次跑在同一 `PackState` 上，`Placed` 迁移序列、`softened` 之和与
      `converged` 与测试内硬编码字面量完全相等；无真实算法、无第三方运行时。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002 / ac-004** 是改名与瘦身的完整性门：四个方法、无 `validate`、无 `Layers`、无别名。
- **ac-003** 直接编码 law P7：阶段不自报裁决，handler 从状态读数。
- **ac-005** 把可重入契约变成可执行证据。
- **ac-006 / ac-007** 是 handler 可见性与文档门（含预算明账）。
- **ac-008** 保证 seam 的测试只验 seam 自己的行为（law § 11）。
- **ac-009 / ac-010** 是数值门与本 spec 的回归场景。
