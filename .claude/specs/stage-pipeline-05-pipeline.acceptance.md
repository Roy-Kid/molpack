---
slug: stage-pipeline-05-pipeline
criteria:
  - id: ac-001
    summary: Pipeline owns the lifecycle body; run is required and has one spelling
    type: code
    pass_when: |
      `src/pipeline/engine.rs` 定义 `pub struct EngineSetup<'a>`、`pub trait StageFactory`
      （含 provided `fn take_handlers(&mut self) -> Vec<Box<dyn Handler>>` 与
      `impl<T: StageFactory + ?Sized> StageFactory for Box<T>`）与
      `pub trait PackEngine: StageFactory + Sized`，其中 `fn run(...)` 是**必需**方法
      （签名以 `;` 结束，无 provided 体），文件 ≤ 300 行；
      `src/pipeline/mod.rs` 定义 `pub struct Pipeline` 与 `new` / `single` / `with_stage`
      并 `impl StageFactory` + `impl PackEngine`，文件 ≤ 400 行；
      `grep -rn 'fn execute' src/pipeline/` 无命中；
      `grep -rn 'fn solver(\|fn prepare(' src/` 无命中。
    status: pending

  - id: ac-002
    summary: One direction only — entry never names pipeline, pipeline never names an algorithm
    type: code
    pass_when: |
      `grep -rn 'pipeline' src/entry/` 无命中；
      `grep -n 'EngineSetup' src/entry/` 无命中（该类型已迁至 `src/pipeline/engine.rs`）；
      `grep -rn 'use crate::\(grow\|gencan\|refine\|initial\)' src/pipeline/` 无命中；
      `grep -n 'use crate::initial' src/entry/mod.rs` 无命中；
      `src/entry/mod.rs` 的模块文档把 `entry` 描述为设置 / 空间 / 结果，不再自称拥有生命周期；
      `src/lib.rs` 在 crate 根重导出 `PackEngine` / `PackSettings` / `PackResult`
      （公开路径 `molpack::PackEngine` 等与本 spec 之前一致）；
      `src/objective.rs`、`src/context/pack_context.rs`、`src/gencan/mod.rs` 行数相对
      本 spec 之前未增。
    status: pending

  - id: ac-003
    summary: Stage-order violations are named before any handler is notified
    type: runtime
    pass_when: |
      `tests/pipeline.rs` 中一个 `requires().placed == Placed::All` 的假阶段排在最前时返回
      `PackError::StageOrder { stage, needs }`（两字段均为 `&'static str`），且一个记录
      `on_start` 的 handler 断言未被调用、无任何阶段运行；该错误的 `Display` 同时点名阶段名、
      缺失的前置条件与补救办法。
    status: pending

  - id: ac-004
    summary: A preset carrying non-default shared settings into with_stage is refused by name
    type: runtime
    pass_when: |
      `Pipeline::new().with_stage(GenCanPack::new().with_seed(7)).run(..)` 返回
      `PackError::PresetSettingsInsidePipeline { stage, knob }`，`knob` 为 `"seed"`，
      `Display` 点名阶段与旋钮并给出补救办法；`GenCanPack::new().with_seed(7).run(..)`
      （走 `Pipeline::single`）**不**触发该错误；
      `PackSettings::first_non_default_knob` 的实现体对 `PackSettings` 做完整解构
      （`grep -n '\.\.' src/entry/mod.rs` 在该函数体内无命中）。
    status: pending

  - id: ac-005
    summary: A preset's handlers are adopted by the pipeline, never dropped
    type: runtime
    pass_when: |
      `tests/pipeline.rs` 断言 `Pipeline::new().with_stage(GenCanPack::new()
      .with_handler(Box::new(counting)))` 运行后该 handler 的 `on_step` 计数 > 0；
      在 `[CbmcGrow.with_handler(c1), GenCanPack.with_handler(c2)]` 下 `c1` 与 `c2` 都收到
      两个阶段的回调（`StepInfo.stage.index` 取到 0 与 1），采纳顺序为阶段顺序；
      `src/pipeline/engine.rs` 的 `take_handlers` rustdoc 与 `src/pipeline/mod.rs` 的
      `with_stage` rustdoc 都写明"handler 被采纳、settings 被具名拒绝"。
    status: pending

  - id: ac-006
    summary: Single-stage pipeline equals the direct run, bitwise
    type: runtime
    pass_when: |
      `tests/pipeline.rs` 断言 `Pipeline::new().with_stage(X).run(t, n)` 与 `X.run(t, n)` 的
      `positions()` 与 `fdist` / `frest` 逐位相同（`to_bits()`），
      X ∈ {`GenCanPack`, `CbmcGrow`, `LatticeGrow`}；其中无 box / cell 声明的 `GenCanPack`
      也在集合内（`Placed::None` 且无种子时前奏不装网格、不 panic）。
    status: pending

  - id: ac-007
    summary: Grow-then-GENCAN in a pipeline equals seeded_from, bitwise
    type: runtime
    pass_when: |
      `tests/pipeline.rs` 断言
      `Pipeline::new().with_stage(CbmcGrow::new(p)).with_stage(GenCanPack::new()).run(t, n)`
      与 `GenCanPack::new().seeded_from(&cbmc_result).run(t, n)`（同 seed、同共享设置）的
      `positions()` / `fdist` / `frest` 逐位相同（`to_bits()`）；
      `src/pipeline/mod.rs` 在每个阶段 `run` 之前调用 `state.invalidate_geometry_cache()`
      （`grep -n 'invalidate_geometry_cache' src/pipeline/mod.rs` 命中该边界调用），
      使两条拼写对称冷启动。
    status: pending

  - id: ac-008
    summary: A second GENCAN stage continues from the first — it never re-initialises
    type: runtime
    pass_when: |
      `tests/pipeline.rs` 断言
      `Pipeline::new().with_stage(GenCanPack::new()).with_stage(GenCanPack::new()).run(t, n_small)`
      与 `GenCanPack::new().seeded_from(&first_result).run(t, n_small)`
      （`first_result = GenCanPack::new().run(t, n_small)`，同 seed、同共享设置，`n_small`
      小到第一阶段不收敛）的 `positions()` / `fdist` / `frest` 逐位相同（`to_bits()`）；
      `grep -rn 'push_off' src/` 在 `GencanSettings` 上无字段命中，push-off 判据由
      `state.placed() == Placed::All` 派生。
      该场景用**声明了周期盒**的输入（两条拼写解析出同一个 cell）；比较只针对
      `positions()` / `fdist` / `frest`，不比较 `frame.simbox`（无盒的管线不设它，`seeded_from`
      继承并设置它）。前奏在两条拼写下都以 `Placed::All` 安装同一个盒、同一个 `radmax`。
    status: pending

  - id: ac-009
    summary: One verdict, one Placements derivation, honest early stop
    type: runtime
    pass_when: |
      `grep -rn 'fdist:' src/entry/ src/pipeline/` 显示 `PackResult.fdist` 只由末阶段结束后的
      状态填充（无第二处赋值、无额外评估调用、无对 `StageOutcome` 的读取）；
      `Placements` 恒从状态的刚体视图槽派生（`src/pipeline/mod.rs` 内无按阶段类型分支的代码）；
      `StepInfo.stage.index` 单调且 `total` 等于阶段数；`on_stage_start` / `on_stage_end`
      各被调用阶段数次；中途 `should_stop` 后续阶段不运行、`converged == false`、键长仍等于
      模板值（1e-9）。
    status: pending

  - id: ac-010
    summary: The Python gate stays green across the trait move
    type: runtime
    pass_when: |
      `python/src/entry.rs` 使用 `use molpack::PackEngine;`（`grep -n 'molpack::entry::PackEngine'
      python/src/` 无命中）；`uv run --directory python --group dev tox -e py` 全绿，
      测试集合与本 spec 之前一致（本 spec 不新增 Python 测试）。
    status: pending

  - id: ac-011
    summary: The Packmol regression and the suite stay green
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` 通过（既有 RED
      `grow_cg_kremer_grest_c_inf` 除外，未被修改或跳过）；
      `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored`
      五例通过；`seeded_run_contract` 与 `free_chain_push_off_deterministic` 断言未改且全绿。
    status: pending

  - id: ac-012
    summary: Regression scenario reproduces the hard-coded pipeline goldens
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- pipeline_regression_single_stage_gencan_golden`
      通过：固定输入下的 `fdist` 与前三个原子坐标与测试内硬编码字面量在 1e-12 内相等；
      注释记录金标捕获自本 spec 之前的构建；无第三方运行时。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** 是架构门：生命周期只有一个体、一个动词，依赖只有一个方向，`EngineSetup`
  与它的生产者 / 消费者同住。
- **ac-003 / ac-004 / ac-005** 是"不静默"门：两类用户错误在跑之前具名报出，新旋钮无法逃过检查，
  预设随身的 handler 被采纳而不是被丢掉。
- **ac-006 / ac-007 / ac-008** 是三道数值门：预设逐位不变；管线拼写与 `seeded_from` 拼写在
  跨算法与同算法两种衔接下都重合（law P7 的"显式串联"两种写法必须一致），后者同时证明第二个
  同算法阶段不会重新 `initial()`。
- **ac-009** 钉住单一裁决、单一 `Placements` 派生与诚实的 handler / early stop。
- **ac-010 / ac-011** 是两条既有门（Python 与 Packmol 回归）。
- **ac-012** 是本 spec 的回归场景。
