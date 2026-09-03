---
title: stage-pipeline-05-pipeline — Pipeline：生命周期体的唯一所在
status: approved
created: 2026-09-02
chain: stage-pipeline（05 of 7）
---

# stage-pipeline-05-pipeline

## Summary

引入 `Pipeline`（`src/pipeline/mod.rs`）：持有共享设置、handler 与一串阶段，是生命周期**体**的唯一
所在（校验 → 建状态 → 阶段衔接检查 → 逐阶段运行 → 裁决 → 组装）。`StageFactory`、`PackEngine`
与 `EngineSetup` 一并迁入 `src/pipeline/engine.rs`——生命周期的所有者持有生命周期的 trait 与它的
入参；`src/entry/` 收缩为设置 / 空间 / 结果（`PackSettings`、`setup.rs`、`result.rs`），**永不提及**
`pipeline`。`PackEngine::run` 成为必需方法，三个预设各以一行实现
（`Pipeline::single(self).run(targets, max_loops)`），`Pipeline` 自己的 `run` 就是生命周期体。
push-off 从 `GencanSettings` 的旗标升为状态词汇（`Placed::All`），因此**任何**前驱阶段之后的
GENCAN 阶段都从既有放置继续而不是重新 `initial()`；进入管线的预设携带的 handler 被采纳而不是被
丢弃。`GenCanPack` / `CbmcGrow` / `LatticeGrow` 的公开 API 与坐标逐位不变；衔接错误在任何 handler
被通知之前具名报出。

## Design

**层次与依赖方向**（评审 🚨 的直接落实）。上一稿让 `entry` 的 provided `run` 去包 `Pipeline`，
同时 `pipeline` 又必须认识 `entry` 的一切——一条闭环。本稿只留一个方向：

    pipeline → { entry(settings/space/result), stage, context, objective, handler, error, target }
    entry    → { context, handler, target, error }        （不含 pipeline）

- `src/pipeline/engine.rs`（new，预算 ≤ 300 行）持 `StageFactory`、`PackEngine`（含全部 provided
  `with_*` builders，逐字自 `src/entry/mod.rs:108-215` 迁入）与 `EngineSetup`（自
  `src/entry/mod.rs:89-98` 迁入）。
- **`EngineSetup` 跟着它的生产者与消费者走**（评审 🟡 的落实）：它由生命周期体填充
  （`src/entry/mod.rs:248-257`，本 spec 搬进 `Pipeline::run`），只被 `StageFactory::stages` 消费；
  生产者与消费者都在 `pipeline/` 之后，把类型留在 `entry/` 就是一件事两个所有者（law § 4），
  而且会让仍叫 `entry` 的模块继续假装拥有"入口"这个概念。`src/entry/mod.rs` 的模块文档改写为
  **entry = 共享设置 + 空间解析 + 结果**（`PackSettings` / `LogSpec` / `setup.rs` / `result.rs`），
  入口本身住在 `gencan/` / `grow/`，生命周期住在 `pipeline/`。
- `src/pipeline/mod.rs`（new，预算 ≤ 400 行）持 `Pipeline` 与生命周期体（逐段自
  `src/entry/mod.rs:216-346` 迁入）。
- `src/entry/mod.rs` 保留 `LogSpec` / `PackSettings` 与 `result` / `setup` 的模块声明，并新增
  `PackSettings::first_non_default_knob`。它对 `pipeline` 零引用——连 `pub use` 都不写（否则
  grep 门会命中）。
- **公开面不变**：`src/lib.rs` 继续在 crate 根重导出 `PackEngine` / `PackSettings` / `PackResult`
  （新增 `StageFactory` / `Pipeline`），用户今天写的 `use molpack::PackEngine;` 一字不改。
  失效的是两条**模块限定**路径 `molpack::entry::PackEngine` 与 `molpack::entry::EngineSetup`
  （后者今天经 `pub mod entry` 可达，本 spec 后为 `molpack::pipeline::EngineSetup`）；仓库内的
  使用者只有三个预设入口与 `python/src/entry.rs:12`，都在本 spec 的 Files 里跟随修改。这是一次
  路径迁移，在此明账（law § 10），不做兼容别名（law § 1）。

**`StageFactory` 与 `PackEngine`**（`src/pipeline/engine.rs`）：

    pub struct EngineSetup<'a> { /* 字段逐字自 entry/mod.rs:89-98 */ }

    pub trait StageFactory {
        fn validate_targets(&self, targets: &[Target]) -> Result<(), PackError> { Ok(()) }
        fn settings(&self) -> &PackSettings;
        /// Surrender the handlers this factory carries. Default: none.
        fn take_handlers(&mut self) -> Vec<Box<dyn Handler>> { Vec::new() }
        fn stages(&mut self, setup: &EngineSetup<'_>) -> Result<Vec<Box<dyn Stage>>, PackError>;
    }
    impl<T: StageFactory + ?Sized> StageFactory for Box<T> { /* 四个方法逐一转发 */ }

    pub trait PackEngine: StageFactory + Sized {
        fn settings_mut(&mut self) -> &mut PackSettings;
        fn handlers_mut(&mut self) -> &mut Vec<Box<dyn Handler>>;
        // …全部 provided `with_*`（原样迁入）…
        fn run(self, targets: &[Target], max_loops: usize) -> Result<PackResult, PackError>;
    }

- `solver()` → `stages()`（返回向量，预设各返回一个阶段）；`prepare()` **删除**——它的两件事分别
  下沉到 `GencanStage::run` 的前奏（见下）。
- **`take_handlers` 是 `StageFactory` 的（评审 🔴 的落实）**。handler 今天挂在 `PackEngine` 上
  （`handlers_mut`，`entry/mod.rs:166-172`），而 `with_stage` 只看得见 `StageFactory`，于是
  `Pipeline::new().with_stage(GenCanPack::new().with_handler(h))` 会让 `h` 静默失灵——恰恰是
  `entry/mod.rs:8-11` 说这套设计要防的那件事。修法：在 `StageFactory` 上加一个 provided
  `take_handlers`，默认交出空集；三个预设各以 `std::mem::take(self.handlers_mut())` 实现
  （`PackEngine: StageFactory`，`handlers_mut` 就在手边）；`Pipeline` 自己覆写它交出自己的
  handler 集，因此嵌套管线也不丢。
- **`impl StageFactory for Box<T>`**：07 的 `IntoStageFactory` 交出的是
  `Box<dyn StageFactory>`，`with_stage` 必须能直接吃下已装箱的工厂。转发实现在此写明，
  不留给实现期发现（law § 7）。
- **`run` 是必需方法**（评审 🟡 (d) 的落实）：没有 provided 体，就没有"把自己再包一层"的隐患，
  也就不需要覆写与解释覆写。三个预设各一行：
  `fn run(self, t, n) -> … { Pipeline::single(self).run(t, n) }`。
  `Pipeline` 直接 `impl PackEngine`，它的 `run` **就是**生命周期体。
- **一个动词**：`Pipeline` 没有固有 `execute`。同一件事只有 `run` 这一个拼写（law § 7）。

**`Pipeline`**（`src/pipeline/mod.rs`）：

- `Pipeline::new()` + `with_stage(impl StageFactory + 'static)`：**组合**拼写。
  **命名**：crate 里约 30 个消费型 builder 无一例外是 `with_*`（`entry/mod.rs:169 with_handler`、
  `target.rs:387 with_restraint`），父 spec 的 `.stage(...)` 未给理由，本 spec 统一为 `with_stage`；
  签名形状照 `Target::with_restraint`（调用点 `impl Trait + 'static`，内部装箱）。
- **`with_stage` 对两类随身之物的两种处置，都不是静默丢弃（law § 8 / § 10）**：
  - **handler → 采纳**。`self.handlers.extend(stage.take_handlers())`，按**阶段顺序**追加进管线
    的 handler 集；rustdoc 写明这一点。理由：handler 是观察者，跨阶段观察本来就是它的语义，
    采纳不会改变任何一把尺子。
  - **非默认共享 `PackSettings` → 具名拒绝**。`PackSettings` 是**一把尺子**（tolerance /
    precision / seed / cell …），两个阶段各带一把，共享 objective 就不再唯一——所以按名报错
    `PackError::PresetSettingsInsidePipeline`，而不是挑一个赢家。
- `Pipeline::single(engine: impl PackEngine + 'static)`：**预设**拼写。它采纳该 engine 的
  `PackSettings`，并同样经 `take_handlers()` 采纳其 handler 向量（与 `with_stage` 一个拼写，
  不是两套排空逻辑），再把 engine 作为唯一阶段来源装箱。它与 `with_stage` 的差别只在设置：
  `single` 采纳设置，`with_stage` 拒绝非默认设置。
- `Pipeline` 亦 `impl StageFactory`（展平自己的阶段、覆写 `take_handlers`），因此可嵌套；
  `settings()` 返回自己的设置。
- 生命周期体逐段搬迁自 `entry/mod.rs:216-346`：① 空目标 / 空分子校验 + 每个 factory 的
  `validate_targets` → 全局约束广播 + 空间解析 → ② `build_context` + `PackState` 组装 →
  **阶段衔接检查** → 逐阶段 `run`（每阶段前失效几何缓存，前后调 `on_stage_start` /
  `on_stage_end`，返回后按 `guarantees().placed` 推进状态的 `placed`）→ 裁决 → ⑤ 组装。
- **阶段衔接检查**：把 `Placed` 沿链推进——初始 `Placed::None`，每个阶段跑完按其
  `guarantees().placed` 推进；某阶段 `requires().placed == All` 而当前是 `None` →
  `PackError::StageOrder`。检查在 `stages()` 解析之后、**任何 handler 被通知之前、任何阶段运行
  之前**返回。（不把 `requires` / `guarantees` 复制到 `StageFactory` 上——那是同一事实两个家。）
- **错误形状**（评审 🟡 的落实，两者的载荷都是**纯字符串**，`error` 不获得指向 `context` /
  `stage` 的新边）：

      StageOrder { stage: &'static str, needs: &'static str }
      PresetSettingsInsidePipeline { stage: &'static str, knob: &'static str }

  `needs` 是渲染后的前置条件名（如 `"placed: all"`）。`Display` 点名冒犯者**与**补救办法——
  `SeedMismatch`（`error.rs:110`）与 `UnknownMass`（`:105`）是最近的模型。
- **预设带非默认共享设置进 `with_stage`** 的检测经
  `PackSettings::first_non_default_knob(&self) -> Option<&'static str>`，实现方式是**完整解构**：

      let PackSettings { tolerance, precision, discale, seed, parallel_eval,
                         short_tolerance, periodic_box, density, cell, log,
                         global_restraints } = self;   // 没有 `..`

  新增一个字段会**编译失败**，而不是静默逃过检查（评审 🟡 的落实；`PackSettings` 含
  `Vec<Arc<dyn AtomRestraint>>` 无法 `PartialEq`，故逐项判定并回报第一个越界旋钮名）。

**`GencanStage::run` 的前奏与 push-off 的状态化**（评审两条 🔴 + 两条 🟡 的落实）。`prepare()` 消失后，
它做的两件事有了明确的落点，且判据全部用状态词汇表达：

    GencanStage::run(state, targets, budget, handlers):
      ① 网格 / simbox：**仅当本阶段将从既有放置继续时**安装——即
         `state.placed() == Placed::All`，或本预设携带种子（`seeded_from`，它在 ② 把状态置为
         `Placed::All`）。盒子取 `setup.cell`（已解析时）否则取上下文当前的 `ctx.simbox`；
         `radmax` 与今天 `src/gencan/entry.rs:185` 一致取 `max(sys.radius)`
         （`install_simbox_and_grid(ctx, cell, radmax, self.discale, ntotat_free)`，
          逐字自 `:181-192`）。`Placed::None` 且无种子时**不装**：`initial()` 自己拥有盒子与
         网格（无声明时用 `sidemax` 合成回退盒，`src/initial.rs:536-549`；有声明时 `:549` 重装）。
         这样两条拼写在同一状态下装的是同一个盒、同一个 `radmax`，
         `[GenCanPack, GenCanPack] ≡ seeded_from(first)` 在无盒场景也成立。今天 `:181-184` 的
         `.expect("seeded_from installed the cell declaration")` 不随行。
         （二次安装是幂等的：`resize_cell_arrays` 会清 `latomfix` / `fixed_cells`，
          `src/context/pack_context.rs:408-421`。）
         **既有债务，在此明账（law § 10）**：push-off 路径的 `radmax = max(sys.radius)`
         （`gencan/entry.rs:185`、`grow/entry.rs:138`）与 `initial()` 的
         `radmax = 2·max(radius_ini)`（`src/initial.rs:456-461`）是同一事实的两种推导，
         前者的 ±1 模板覆盖约 `1.01·discale·R`，而配对截断是 `2·discale·R`。本 spec 只保证
         两条拼写用**同一种**推导（因而逐位一致），不修正推导本身——修正归 `/mol:debug`
         （债务记录 D-02，`.claude/notes/notes.md`）。
      ② **种子注入点**（本 spec 唯一的一处）：若本预设携带种子（`GenCanPack::seeded_from`），
         `*view = RigidView::install_seed(seed.rigid.as_slice(), &seed.coor, ctx)`，
         随即 `state.set_placed(Placed::All)`。顺序是先 ①（装网格）后 ②（拷坐标），
         与今天 `prepare()` 内的顺序一致——逐位等价依赖这个顺序。
      ③ push-off 判据：`let push_off = state.placed() == Placed::All;`
         为真 → 从既有放置继续：跳过 `initial()`，`RigidView::write_xcart` 物化 xcart，
         movebad 关闭（今天 `src/gencan/solver.rs:156-162` 的语义逐字保留）；
         为假 → `initial()`（今天 `:143-155`）。

- **`GencanSettings.push_off` 字段在本 spec 删除**（`src/gencan/solver.rs:51,66,83,98,105`），
  `src/gencan/entry.rs:213` 那行推导随之消失。`:143` 的分支与 `:215` 传给 `run_phase` 的
  `self.push_off || !self.settings.perturb` 都改读 ③ 的派生局部量。事实只有一个家
  （`PackState.placed`），没有第二份可变表示（law § 9）。02 把该字段留在原处不动，本 spec 是它
  唯一的删除点。
- `GencanStage::requires()` 保持 `Placed::None`（它能从零开始），`guarantees()` 是 `Placed::All`。
  管线在每个阶段 `run` 返回后按 `guarantees().placed` 推进状态。因此 `[X, GenCanPack]` 中的 GENCAN
  阶段**一定**看见 `Placed::All`——`X` 是生长、是另一个 GENCAN、还是 06 的组合子都一样，没有哪条
  路径会把前一阶段的成果 `initial()` 掉。生长阶段 `run` 末尾的 `capture_from_xcart`（02）保证槽里
  的视图有效，所以 `Placed::All` 同时就是"槽里视图有效"的说法。

**阶段边界的几何缓存**（评审 🟡 的落实）。管线在**每个**阶段 `run` **之前**调
`state.invalidate_geometry_cache()`（03 的转发器，其唯一调用方就在这里）。理由：一条管线复用同一个
`PackContext`，前一阶段末尾的 `evaluate_unscaled` 会把几何缓存留成**热**的（键 = x + comptype +
init1 + `geometry_key`，`src/context/work_buffers.rs:75-88`），而 `seeded_from` 拼写是从**冷**缓存
出发的；两个前奏用相同的 `radmax` / `discale` 装网格，于是命中与否会让
`src/objective.rs:140-146` / `:180-187` 走**不同的求和路径**（命中会跳过 `resetcells()` +
`expand_molecules`，改经 `accumulate_constraint_values_from_xcart` 累加约束），`frest` 的逐位相等
就不再有保证。在每个边界失效使两条拼写对称地冷启动。第一个阶段之前上下文本就是新建的，这一次调用是
空操作，因此单阶段的逐位等价（ac-006）不受影响。

**两组逐位等价**（本 spec 的核心数值断言）。前提逐条列出并由测试锁定：两侧 `coor` / `x` 同源
（02 的同一写回契约与同一 `install_seed`）；两侧网格由 GENCAN 前奏用相同参数安装；两侧
`scale` / `scale2` 都在默认值（03 的对称存还保证前一阶段结束时已还原）；两侧 push-off 为真故
movebad 关闭、RNG 由同一 `settings.seed` 播种；`cell` 由同一份共享设置解析；两侧几何缓存都冷。

    Pipeline::new().with_stage(CbmcGrow::new(p)).with_stage(GenCanPack::new()).run(t, n)
      ≡（逐位：positions / fdist / frest）
    GenCanPack::new().seeded_from(&cbmc_result).run(t, n)          // 同 seed、同 settings

    Pipeline::new().with_stage(GenCanPack::new()).with_stage(GenCanPack::new()).run(t, n_small)
      ≡（逐位：positions / fdist / frest）
    GenCanPack::new().seeded_from(&first_result).run(t, n_small)
      // first_result = GenCanPack::new().run(t, n_small)；n_small 刻意取小使第一阶段不收敛
      // （v1 一个 Budget 给全部阶段，两阶段同预算）

第二组证明"**同算法**的第二个阶段从既有放置继续，绝不重新 `initial()`"——这是 ③ 的可观察形式，也是
06 的 `Repeat` 赖以成立的前提。它的第一阶段与 `first_result` 的逐位相同由单阶段等价（ac-006）保证。

**`Placements` 的派生**（评审 🟡 的落实，只用状态词汇，不问算法身份）：写回契约是——**每个阶段
结束时必须让状态里的刚体视图有效**（GENCAN 写回它的 `x`；生长阶段 `capture_from_xcart`，02 已
落地）。因此管线**恒**从 `state` 的槽里取 `Placements`，没有任何分支。rustdoc 把 `seeded_from`
的逐位连续性限定为"最后一个阶段维护了刚体视图的运行"——今天三个预设与本链所有阶段都满足。

**裁决**：`PackResult.fdist` / `.frest` 取自**最后一个阶段运行结束后的状态**（`state.ctx().fdist` /
`.frest`），即今天 `entry/mod.rs:341-342` 的语义，**不额外做一次评估**。这在 03 之后是安全的：
每个阶段的最后一个动作都是那一个未缩放裁决原语。管线永不读 `StageOutcome`（04 已把裁决字段从
它上面拿掉）。`converged = 末阶段 converged && fdist < precision && frest < precision`；
`softened` = 各阶段之和。

**handler**：`on_start` / `on_finish` 每 run 一次（不是每阶段一次）；`on_stage_start` /
`on_stage_end` 每阶段一次；`StepInfo.stage.index` 单调、`total` 等于阶段数；任一阶段内
`should_stop` 生效即终止后续阶段，`converged = false`，坐标仍经该阶段的 abort 契约保持成键。
被采纳进来的预设 handler 与管线自己的 handler 在同一个集合里，看得见**全部**阶段。

### Reuse decision

- `entry/mod.rs:108-215 PackEngine` + `:169 with_handler` + `handlers_mut()` — **generalize**：整体迁入 `src/pipeline/engine.rs`；`Vec<Box<dyn _>>` 的拥有 + 一次性排空（`:268 std::mem::take`）升为 `StageFactory::take_handlers` 的预设实现，`single` 与 `with_stage` 共用这一个拼写。
- `entry/mod.rs:89-98 EngineSetup` — **generalize**：随它的生产者（生命周期体）与消费者（`StageFactory::stages`）迁入 `src/pipeline/engine.rs`，`entry/` 不再持有它。
- `entry/mod.rs:216-346` 生命周期体 — **generalize**：整体搬入 `src/pipeline/mod.rs`，成为 `Pipeline` 的 `PackEngine::run`。
- `target.rs:387 with_restraint` — **pattern**：`with_stage(impl StageFactory + 'static)` 的调用点签名。
- `error.rs:5-56 PackError` — **pattern**：两个新变体具名字段 struct-like（纯 `&'static str` 载荷）+ 点名补救的 `Display`；`DensityConflictsWithBox`（`:48`）是同形状但无数据的单元变体，故不照抄其形。
- `entry/mod.rs:341-342 fdist: sys.fdist` — **reuse**：裁决口径保持"读状态"，不新增评估；这正是逐位不变的前提。
- `gencan/entry.rs:170-192 prepare` 的网格 / simbox 安装 — **reuse**：下沉为 `GencanStage::run` 前奏的 ①，**条件化**为 `state.placed() == Placed::All` 或携带种子（盒子取 `setup.cell` 否则 `ctx.simbox`），`.expect(...)` 不随行（无盒回退仍归 `initial()`）；`radmax` 两种推导的分歧记为债务 D-02。
- `context/rigid_view.rs::{install_seed, capture_from_xcart, write_xcart}`（02）— **reuse**：种子注入、写回契约与 xcart 物化的唯一来源，管线不重新推导。
- `gencan/solver.rs:51 GencanSettings.push_off` — **generalize**：字段在本 spec 删除，判据升为状态词汇 `state.placed() == Placed::All`（02 把它留在原处不动）。
- `context/pack_state.rs::{evaluate_unscaled, invalidate_geometry_cache}`（03）— **reuse**：裁决原语不重复；缓存失效在每个阶段边界调用，这是它的唯一调用方。
- `src/entry/result.rs::Placements` — **reuse**：形状由 02 定，管线只从状态槽填充它。

## Files to create or modify

- `src/pipeline/mod.rs` (new)
- `src/pipeline/engine.rs` (new)
- `src/entry/mod.rs`
- `src/entry/result.rs`
- `src/error.rs`
- `src/gencan/entry.rs`
- `src/gencan/solver.rs`
- `src/grow/entry.rs`
- `src/grow/driver.rs`
- `src/grow/lattice/entry.rs`
- `src/grow/lattice/mod.rs`
- `src/lib.rs`
- `python/src/entry.rs`
- `tests/pipeline.rs` (new)

## Tasks

- [ ] Write failing tests in `tests/pipeline.rs`：单阶段等价 ×3（`Pipeline::new().with_stage(X).run(..)` 与 `X.run(..)` 的 `positions()` / `fdist` / `frest` 逐位相同）、`[CbmcGrow, GenCanPack]` ≡ `GenCanPack::seeded_from(cbmc_result)` 逐位、`[GenCanPack, GenCanPack]`（小 `max_loops`）≡ `GenCanPack::seeded_from(first_result)` 逐位、`Pipeline::new().with_stage(GenCanPack::new().with_handler(counting))` 下该 handler 收到 `on_step` 回调、`StageOrder` 在任何 handler 被通知前触发、`PresetSettingsInsidePipeline` 点名旋钮、两阶段管线的 `softened` 求和、`StepInfo.stage` 单调与 `on_stage_start/end` 计数、中途 `should_stop` 保持成键几何且 `converged == false`；确认 RED
- [ ] Add `src/pipeline/engine.rs`：`EngineSetup`（自 `src/entry/mod.rs:89-98` 迁入）、`StageFactory`（含 provided `take_handlers` 与 `impl<T: StageFactory + ?Sized> StageFactory for Box<T>`）与 `PackEngine`（`run` 为**必需**方法，provided `with_*` 自 `:108-215` 迁入），删除 `solver()` 与 `prepare()`；`src/entry/mod.rs` 只留 `LogSpec` / `PackSettings` 与模块声明并改写模块文档为"设置 + 空间 + 结果"；`src/lib.rs` 在 crate 根重导出 `PackEngine` / `StageFactory` / `Pipeline`
- [ ] Add `src/pipeline/mod.rs`：`Pipeline` + `new` / `single` / `with_stage`（`with_stage` 经 `take_handlers` **采纳**阶段 handler，按阶段顺序追加），`impl StageFactory`（覆写 `take_handlers`）与 `impl PackEngine`（其 `run` 即生命周期体，自 `src/entry/mod.rs:216-346` 迁入，含衔接检查、每阶段前的 `invalidate_geometry_cache`、`guarantees().placed` 推进与两个 stage 钩子），并补模块文档（生命周期五段、衔接检查时机、handler 采纳规则、缓存边界、裁决取自状态、`seeded_from` 连续性的适用范围、写回契约）
- [ ] Implement `StageFactory` on the three presets and their one-line `run`（`src/gencan/entry.rs`、`src/grow/entry.rs`、`src/grow/lattice/entry.rs`：`take_handlers` 用 `std::mem::take(self.handlers_mut())`，`run` 用 `Pipeline::single(self).run(t, n)`），并把网格 / simbox 安装下沉为各 `run` 前奏——GENCAN 侧**仅当 `state.placed() == Placed::All` 或携带种子**时安装，盒子取 `setup.cell` 否则 `ctx.simbox`（`src/gencan/solver.rs`），生长侧照旧（`src/grow/driver.rs`、`src/grow/lattice/mod.rs`）
- [ ] Add `PackError::StageOrder { stage, needs }` and `PackError::PresetSettingsInsidePipeline { stage, knob }` to `src/error.rs`（载荷均为 `&'static str`，`Display` 点名阶段、缺失前置条件 / 越界旋钮与补救办法）与 `PackSettings::first_non_default_knob` in `src/entry/mod.rs`（**完整解构，不写 `..`**）
- [ ] Move push-off into state terms in `src/gencan/solver.rs`：删除 `GencanSettings.push_off` 字段与 `src/gencan/entry.rs:213` 的推导，`run` 前奏做种子注入（`*view = RigidView::install_seed(..)` + `set_placed(Placed::All)`）并派生 `push_off = state.placed() == Placed::All`，`:143` 的分支与 `:215` 传给 `run_phase` 的参数改读它；管线在 `src/entry/result.rs` 恒从状态槽派生 `Placements`（无算法身份分支）与 `softened` 求和
- [ ] Follow the moved trait in `python/src/entry.rs`：`use molpack::entry::PackEngine` 改为 `use molpack::PackEngine`（仅 import 跟随，无绑定面变更）
- [ ] Add regression scenario `pipeline_regression_single_stage_gencan_golden` to `tests/pipeline.rs`（硬编码 `fdist` 与前三个原子坐标金标，容差 1e-12）
- [ ] Verify the release gate and the Python gate：`cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored` 五例、`uv run --directory python --group dev tox -e py`
- [ ] Run full check + test suite

## Testing strategy

- 归属：管线的生命周期契约归 `tests/pipeline.rs`（law § 11）。单测门：
  `cargo test -p molcrafts-molpack --lib --tests -- pipeline`。
- Happy path：单阶段等价 ×3；两阶段管线（`CbmcGrow` → `GenCanPack`）跑通并 `softened` 求和。
- **两条 push-off 数值门（本 spec 的核心断言）**：
  (1) `Pipeline::new().with_stage(CbmcGrow::new(p)).with_stage(GenCanPack::new()).run(t, n)` 与
  `GenCanPack::new().seeded_from(&cbmc_result).run(t, n)` 的 `positions()` / `fdist` / `frest`
  逐位相同（`to_bits()`）；
  (2) `Pipeline::new().with_stage(GenCanPack::new()).with_stage(GenCanPack::new()).run(t, n_small)`
  与 `GenCanPack::new().seeded_from(&first_result).run(t, n_small)` 同样逐位相同——`n_small` 取小
  使第一阶段不收敛，从而"第二个 GENCAN 阶段确实接着跑而不是重来"是可观察的。两条都依赖阶段边界的
  几何缓存失效使双方冷启动（见 Design）。
- **handler 采纳门**：`Pipeline::new().with_stage(GenCanPack::new().with_handler(counting))` 下
  计数 handler 的 `on_step` 计数 > 0；两阶段时它看得见两个阶段（`stage.index` 出现 0 与 1）。
- Edge cases：空阶段列表 → 具名错误；`StageOrder`（把一个 `requires: Placed::All` 的测试用假阶段
  排在最前）；预设带 `with_seed` 进 `with_stage` → `PresetSettingsInsidePipeline { knob: "seed" }`；
  预设自己的 `run`（走 `Pipeline::single`）**不**触发该错误；**无 box / cell 声明**的单阶段 GENCAN
  管线正常跑完（前奏不装网格、不 panic，回退盒仍由 `initial()` 合成）；中途 `should_stop`（用
  `EarlyStopHandler`）→ 后续阶段不运行、`converged == false`、键长仍等于模板值（复用
  `grow_abort_keeps_bonded_geometry` 的判据，1e-9）。
- Regression scenario：`pipeline_regression_single_stage_gencan_golden`，硬编码金标；对应
  `type: runtime` 验收项。
- 既有逐位门：`tests/gencan.rs`、`tests/packer.rs`、`tests/grow.rs` 确定性组、`examples_batch` 五例。
- Python 门：`uv run --directory python --group dev tox -e py`（law P4）——本 spec 只改一行 import，
  该门必须与本 spec 之前同样全绿。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- `Repeat` / `Guarded` / `Invariant`（06）——本 spec 只做线性序列。
- 每阶段独立预算：v1 一个 `Budget` 给全部阶段（`max_loops` 语义由各阶段自述）；若实测饿死或超支，
  再加 `Pipeline::with_stage_budget`（加法扩展）。
- 并行阶段、DAG。
- `.inp` 语法：script 路径按构造仍是单阶段 GENCAN，不新增关键字（law P2）。
- Python 的 `Pipeline` 镜像、`.pyi` / `_protocols` / 文档（07）；本 spec 在 `python/` 只改一行 import。
- 为 `molpack::entry::{PackEngine, EngineSetup}` 保留兼容别名——不做（law § 1），路径迁移已明账。
- `initial ↔ gencan` 的既有环——记录不动。
