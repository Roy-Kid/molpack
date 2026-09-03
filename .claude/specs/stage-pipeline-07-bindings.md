---
title: stage-pipeline-07-bindings — Python Pipeline 镜像与文档
status: approved
created: 2026-09-02
chain: stage-pipeline（07 of 7，收尾）
---

# stage-pipeline-07-bindings

## Summary

把 `Pipeline` 镜像到 Python wheel：`Pipeline([stage, ...])` / `.with_stage(x)` 接受入口对象，
共享 `with_*` 走既有的 `entry_pymethods!`，衔接错误映射为 `ValueError` 且消息含阶段名；`StepInfo`
向 Python 暴露 `stage`（`index` / `total` / `name`）。阶段对象的识别不再是 `PyPipeline` 里的硬编码
下转列表——每个入口 pyclass **自己**提供转换，`PyPipeline` 只调用它，`TypeError` 文案由同一份注册表
生成；阶段对象随身的 Python handler 也随转换被管线采纳，与 05 的策略同形。同步 `.pyi` 与
`_protocols.py`，更新文档站的五个页面与仓库知识文件，重跑 `/mol:map`。

## Design

**层次**：Python 绑定层（`python/src/`、`python/python/molpack/`）与文档。core 不因绑定而改动
（law § 6：核心语义不由绑定定义）——本 spec 不新增任何 Rust core 符号。

**PyO3**（`python/src/entry.rs`）：

- `PyPipeline` 沿用 `SharedKnobs`（`:28`）+ `clone_ref`（`:44`）+ `py_handlers`（`:39`）+
  `entry_pymethods!`（`:150`）这一整套——共享 `with_*` 只绑定一次，与 `PyGenCanPack`（`:233`）/
  `PyCbmcGrow`（`:371`）/ `PyLatticeGrow`（`:521`）完全同形。
- **一个提取点，不是一张下转清单**（评审 🟡 的落实）。上一稿让 `PyPipeline` 对
  `{PyGenCanPack, PyCbmcGrow, PyLatticeGrow}` 逐个 `downcast`，并把同一张名单在 `TypeError`
  文案里再抄一遍；第四个入口（blueprint 已预告的 `refine/`）会强制修改一个无关模块**加上**它的
  错误字符串（law § 4）。本 spec 改成：

  1. 一个 Rust 侧的转换 trait `IntoStageFactory { fn into_stage_factory(&self, py: Python<'_>)
     -> PyResult<Box<dyn molpack::StageFactory>>; }`，**由 `entry_pymethods!` 在每个入口
     pyclass 处一并 stamp**——转换体与它的入口住在一起。`Box<dyn StageFactory>` 能直接进
     `Pipeline::with_stage`，因为 05 已写下 `impl<T: StageFactory + ?Sized> StageFactory for Box<T>`。
  2. 一处 `stage_entry_registry!(PyGenCanPack => "GenCanPack", PyCbmcGrow => "CbmcGrow",
     PyLatticeGrow => "LatticeGrow");` 宏调用，从这**唯一一张表**同时生成 (a) `PyPipeline`
     使用的按序尝试转换的派发函数，与 (b) `TypeError` 文案里"支持的入口"列表。
  3. `PyPipeline` 自身不含任何具体入口类型名（grep 可查）。新增第四个入口 = 在它自己的 pyclass
     旁加一行注册，`PyPipeline` 与错误文案零改动。
- **handler 的采纳与 05 同形**（评审 🟡 的落实）：转换在交出 `Box<dyn StageFactory>` 之前，把该
  Python 阶段对象的 `py_handlers`（`:39`）逐个包成 `PyHandlerWrapper` 挂到它上面，于是
  `StageFactory::take_handlers` 会把它们交给管线、被管线采纳进自己的 handler 集。
  `Pipeline([GenCanPack().with_handler(cb)])` 因此**会**调用 `cb`——不静默丢弃（law § 8）。
  Python 侧不另立策略，只镜像 05：handler 采纳，非默认共享设置具名拒绝。
- 阶段对象按值收取（`Vec<PyObject>` → 经上述派发转换），与既有的 `targets: Vec<PyTarget>`
  按值收取一致。
- `PackError::StageOrder` / `PresetSettingsInsidePipeline` / `InvariantViolated` 经既有的错误映射
  变成 `ValueError`，消息保留阶段名（Rust 侧 `Display` 已点名，绑定层不重写文案——一个事实一个家）。
- `Stage` / `Invariant` / `Layers` / `Topology` / `PackState` / `RigidView` 保持 **Rust-only**，
  与今天 `Solver` 的定位一致：Python 通过挑入口来挑算法，自定义阶段是 Rust 级扩展点
  （v1 的组合子也不暴露）。
- **CLI**：无改动（law P2：配置面只有 `.inp`，不为管线加关键字）。

**`StepInfo.stage`**（`python/src/handler.rs`）：只读属性 `index` / `total` / `name`，形状照今天
`phase` 的暴露方式。

**类型与协议**：`python/python/molpack/molpack.pyi` 加 `Pipeline` 与 `StepInfo.stage`；
`python/python/molpack/_protocols.py` 同步；`python/python/molpack/__init__.py` 导出。

**文档**：`docs/architecture.md`（模块图 + 依赖方向，含 `src/pipeline/` 与 `entry` 的单向关系；
`:37` 的 `evaluate_unscaled` 从 `phases.rs` 移到 `src/context/pack_state.rs`）、`docs/extending.md`
（`:460` 讲半径互换那段仍把 `evaluate_unscaled` 写在 `phases.rs`，随 03 的重新安家更新，并说明它
现在是 crate 内部原语）、`docs/rust/index.md`（入口 = 单阶段预设，`Pipeline` 是组合）、
`docs/python/guide/packer.md` 新增「Composing stages」一节、`docs/python/api-reference.md` 加
`Pipeline` 与 `StepInfo.stage`。
仓库知识：`CLAUDE.md` 的架构表（`src/solver.rs` → `src/stage.rs`，新增 `src/topology.rs` /
`src/context/pack_state.rs` / `src/context/rigid_view.rs` / `src/pipeline/` / `src/invariant.rs`）、
`.claude/notes/architecture.md`（重跑 `/mol:map`）、`.claude/notes/conventions.md` 的测试一节补上
本链新增的 `tests/*.rs`（`topology` / `context_rigid_view` / `context_pack_state` / `stage` /
`pipeline` / `invariant`）。`.claude/notes/notes.md` 已有 2026-09-02 条目记录「Pipeline 现已 earn，
supersede engine-entry-split 的对应 Non-goal」，落地时核对措辞即可。

### Reuse decision

- `python/src/entry.rs:28 SharedKnobs` / `:44 clone_ref` / `:150 entry_pymethods!` / `:39 py_handlers` — **reuse**：`PyPipeline` 直接接入，不新建平行的共享旋钮机制；转换方法由同一个宏 stamp，`py_handlers` 经转换交给 05 的 `take_handlers`。
- `python/src/handler.rs::PyHandlerWrapper` — **reuse**：阶段对象的 Python 回调照今天的包法包装，绑定层不新造 handler 通道。
- `PyGenCanPack` / `PyCbmcGrow` / `PyLatticeGrow` 的 `targets: Vec<PyTarget>` 按值 marshalling — **pattern**：`Pipeline([stage, …])` 照同一形状。
- `python/src/handler.rs` 的 `StepInfo.phase` 暴露方式 — **pattern**：`stage` 照抄。
- `pipeline/engine.rs::{StageFactory::take_handlers, impl StageFactory for Box<T>}`（05）— **reuse**：handler 采纳与装箱工厂的接受都由 05 提供，Python 侧不复制策略。
- Rust 侧 `PackError` 的 `Display`（05 / 06 已点名阶段、旋钮、不变量）— **reuse**：绑定层只转类型不改文案。
- `python/src/helpers.rs::pack_error_to_pyerr` — **reuse**：三个新变体走既有映射，不新增映射机制。
- 本 spec 不新增任何 Rust core 符号；其余候选项由 01–06 裁决。

## Files to create or modify

- `python/src/entry.rs`
- `python/src/handler.rs`
- `python/src/lib.rs`
- `python/python/molpack/molpack.pyi`
- `python/python/molpack/_protocols.py`
- `python/python/molpack/__init__.py`
- `python/tests/test_pipeline.py` (new)
- `docs/architecture.md`
- `docs/extending.md`
- `docs/rust/index.md`
- `docs/python/guide/packer.md`
- `docs/python/api-reference.md`
- `CLAUDE.md`
- `.claude/notes/architecture.md`
- `.claude/notes/conventions.md`

## Tasks

- [ ] Write failing tests in `python/tests/test_pipeline.py`：两阶段管线（`CbmcGrow` → `GenCanPack`）跑通并返回 `PackResult`；单阶段管线与直接 `run` 的 `positions()` 一致；`Pipeline([GenCanPack().with_handler(cb)])` 下 `cb` 收到 `on_step` 回调（计数 handler，断言计数 > 0）；`StepInfo.stage` 的三字段可见且 `index` 单调；衔接错误与"预设带非默认共享设置进管线"错误映射为 `ValueError` 且消息含阶段名（后者含旋钮名）；未知阶段对象 → `TypeError` 且消息列出注册表里的入口；确认 RED
- [ ] Add the per-entry conversion to `python/src/entry.rs`：`IntoStageFactory` 由 `entry_pymethods!` 在每个入口 pyclass 处 stamp（转换前把该对象的 `py_handlers` 包成 `PyHandlerWrapper` 挂到工厂上，交由 05 的 `take_handlers` 采纳），并加入唯一的 `stage_entry_registry!` 调用（同时生成派发与 `TypeError` 文案）
- [ ] Add the `PyPipeline` pyclass to `python/src/entry.rs`（复用 `SharedKnobs` / `entry_pymethods!` / `py_handlers`，阶段对象按值收取并经注册表派发；类体内不出现任何具体入口类型名）
- [ ] Expose `StepInfo.stage` in `python/src/handler.rs` 并在 `python/src/lib.rs` 注册新类
- [ ] Sync `python/python/molpack/molpack.pyi`, `_protocols.py` and `__init__.py` with `Pipeline` and `StepInfo.stage`
- [ ] Update the docs site：`docs/architecture.md`（模块图 + `:37` 的 `evaluate_unscaled` 位置）、`docs/extending.md`（`:460` 的 `evaluate_unscaled` 随 03 重新安家到 `src/context/pack_state.rs`）、`docs/rust/index.md`、`docs/python/guide/packer.md`（新节「Composing stages」）、`docs/python/api-reference.md`
- [ ] Update repository knowledge：`CLAUDE.md` 架构表、`.claude/notes/architecture.md`（重跑 `/mol:map`）、`.claude/notes/conventions.md` 的测试一节
- [ ] Add regression scenario `test_pipeline_regression_two_stage_golden` to `python/tests/test_pipeline.py`（硬编码 `fdist` 与前三个原子坐标金标，容差 1e-12）
- [ ] Run full check + test suite

## Testing strategy

- 归属：绑定层的契约归 `python/tests/test_pipeline.py`（与既有 `python/tests/test_packer.py` 同形）。
  门：`uv run --directory python --group dev tox -e py`（law P4：非可编辑、隔离，`maturin develop`
  不是门）。
- Happy path：两阶段管线跑通；单阶段管线 ≡ 直接 `run`；`with_seed` / `with_tolerance` 等共享旋钮
  在 `Pipeline` 上可用。
- Edge cases：空阶段列表 → `ValueError`；预设带非默认共享设置进管线 → `ValueError` 且消息含阶段名
  与旋钮名；阶段对象自带的 `with_handler` 回调在管线下**仍被调用**（05 的采纳策略的 Python 面）；
  未知阶段对象 → `TypeError` 且消息中的入口清单来自注册表；handler 回调在多阶段下
  `stage.index` 单调。
- Regression scenario：`test_pipeline_regression_two_stage_golden`，硬编码金标（不导入也不子进程
  调用任何第三方打包工具）；对应 `type: runtime` 验收项。
- 文档门：`cargo doc -p molcrafts-molpack --no-deps` 零警告；五个 docs 页面按 Doc plan 更新。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- Python 侧的自定义 `Stage` / `Invariant` / 组合子：v1 保持 Rust-only 扩展点。
- `Topology` / `PackState` / `RigidView` / `Layers` 的 Python 暴露。
- CLI 与 `.inp` 语法（law P2）。
- molrs ABI 线相关的任何改动（law P10）——本 spec 不触碰 capsule 名或版本门。
- `docs/packmol_parity.md` 与基准页面。
