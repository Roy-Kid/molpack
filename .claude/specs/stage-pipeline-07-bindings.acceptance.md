---
slug: stage-pipeline-07-bindings
criteria:
  - id: ac-001
    summary: Pipeline is available from Python with the shared knobs
    type: code
    pass_when: |
      `python/src/entry.rs` 定义 `PyPipeline`（`#[pyclass(name = "Pipeline")]`）并通过
      `entry_pymethods!` 绑定共享 `with_*`，构造接受 `Pipeline([stage, ...])` 与
      `.with_stage(x)`；`python/src/lib.rs` 注册该类；
      `grep -n 'Pipeline' python/python/molpack/molpack.pyi python/python/molpack/_protocols.py`
      两文件均命中。
    status: pending

  - id: ac-002
    summary: A new entry registers itself; PyPipeline holds no downcast list
    type: code
    pass_when: |
      `python/src/entry.rs` 中每个入口 pyclass 通过 `entry_pymethods!` 获得
      `IntoStageFactory` 的实现，且存在**恰好一处** `stage_entry_registry!` 调用；
      `PyPipeline` 的实现体内不出现 `PyGenCanPack` / `PyCbmcGrow` / `PyLatticeGrow`
      任一标识符，`TypeError` 文案由该注册表生成而非字面量重抄。
    status: pending

  - id: ac-003
    summary: StepInfo exposes the stage triple
    type: runtime
    pass_when: |
      `python/tests/test_pipeline.py` 中一个 Python handler 在多阶段管线下读到
      `info.stage.index` / `.total` / `.name` 三个字段，`index` 单调且 `total` 等于阶段数。
    status: pending

  - id: ac-004
    summary: Stage errors surface as ValueError naming the stage; stage handlers survive
    type: runtime
    pass_when: |
      衔接顺序错误与"预设带非默认共享设置进管线"在 Python 侧均抛 `ValueError`，消息包含
      阶段名（预设设置错误还包含旋钮名，文案与 Rust 侧 `Display` 一致，绑定层未重写）；
      未知阶段对象抛 `TypeError` 并列出注册表中的入口；
      且 `Pipeline([GenCanPack().with_handler(cb)])` 运行后计数回调 `cb` 的 `on_step`
      计数 > 0（阶段对象的 `py_handlers` 经 `IntoStageFactory` → `take_handlers` 被管线
      采纳，与 05 同形，不静默丢弃）。
    status: pending

  - id: ac-005
    summary: The Python gate passes in isolation
    type: runtime
    pass_when: |
      `uv run --directory python --group dev tox -e py` 全绿，含新增的
      `python/tests/test_pipeline.py`；wheel 的 feature 集合未变（仍不含 `io`）。
    status: pending

  - id: ac-006
    summary: Docs and repository knowledge are in sync
    type: docs
    pass_when: |
      `docs/architecture.md` 的模块图含 `src/stage.rs` / `src/topology.rs` /
      `src/context/pack_state.rs` / `src/context/rigid_view.rs` / `src/pipeline/` /
      `src/invariant.rs` 且不再提 `src/solver.rs`，并画出 `pipeline → entry` 的单向关系，
      `:37` 不再把 `evaluate_unscaled` 列在 `phases.rs`；
      `docs/extending.md`（今天的 `:460`）把 `evaluate_unscaled` 指向
      `src/context/pack_state.rs` 并说明它是 crate 内部原语；
      `docs/rust/index.md` 说明"入口 = 单阶段预设、`Pipeline` 是组合"；
      `docs/python/guide/packer.md` 有「Composing stages」一节；
      `docs/python/api-reference.md` 含 `Pipeline` 与 `StepInfo.stage`；
      `CLAUDE.md` 架构表与 `.claude/notes/architecture.md`、`.claude/notes/conventions.md`
      已同步（后者的测试一节列出本链新增的六个 `tests/*.rs`）；
      `cargo doc -p molcrafts-molpack --no-deps` 零警告。
    status: pending

  - id: ac-007
    summary: Regression scenario reproduces the hard-coded Python goldens
    type: runtime
    pass_when: |
      `uv run --directory python --group dev tox -e py -- -k test_pipeline_regression_two_stage_golden`
      通过：两阶段管线的 `fdist` 与前三个原子坐标与测试内硬编码字面量在 1e-12 内相等；
      测试不 import、不子进程调用任何第三方打包工具。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-003 / ac-004** 是绑定面的完整性与错误诚实度；ac-004 同时镜像 05 的 handler
  采纳策略（回调不会因为进了管线而消失）。
- **ac-002** 是 locality 门：第四个入口只改自己那一行注册，`PyPipeline` 与错误文案零改动。
- **ac-005** 是项目规定的唯一 Python 门（law P4）。
- **ac-006** 收口文档与仓库知识，含 `evaluate_unscaled` 两处位置的跟随与 `/mol:map` 刷新。
- **ac-007** 是本 spec 的回归场景。
