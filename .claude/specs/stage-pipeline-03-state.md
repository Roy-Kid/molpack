---
title: stage-pipeline-03-state — PackState：阶段之间交换的状态
status: approved
created: 2026-09-02
chain: stage-pipeline（03 of 7）
---

# stage-pipeline-03-state

## Summary

引入 `PackState`（`src/context/pack_state.rs`）：**包裹**而非抽取 `PackContext`，再加上链式契约需要
的两样东西——放置形状标记 `Placed` 与刚体视图槽 `rigid`。同时把"在 scale = 1.0 的未缩放半径下评估
一次共享 objective"这件今天有两套写法、三处副本的事收进单一实现，三处调用点全部改接，并把
`scale` / `scale2` 的存取改成与 `radius` 同样对称的存—还。GENCAN 与两条生长路径的 `fdist` /
`frest` 逐位不变。本 spec 内 `PackState` / `Placed` 一律 `pub(crate)`。

## Design

**层次**：共享状态底座（`src/context/`）。`gencan/phases.rs` 与两个生长驱动只是被搬迁原语的调用点。

**为什么文件叫 `pack_state.rs`**：`src/context/state.rs` 已被 `RuntimeState`（48 行的**借用**只读遥测
视图，`pack_context.rs:390,396` 发放）占用；`PackState` 是**拥有** `PackContext` 的运行期状态，两者
生命周期故事与职责都不同，不是同一概念的别名。

**类型**（预算 ≤ 300 行）：

    pub(crate) enum Placed { None, All }        // Debug + Clone + Copy + PartialEq + Eq
    pub(crate) struct PackState {               // Debug
        ctx: PackContext,
        placed: Placed,
        rigid: RigidView,
    }

- **可见性**（评审 🟡 (a) 的落实）：`PackState` 与 `Placed` 在本 spec 内是 `pub(crate)`，
  `src/lib.rs` 一行不加。它们在 04 随 `Stage::run(&mut PackState, …)` 的签名成为公开面时才升为
  `pub`。已发布的 `pub` 符号是回退成本最高的部分；在本链尚未落到 seam 之前不付这笔成本。
- **包裹而非抽取**：`PackContext` 的字段一个不动（`src/objective.rs` 的读集不变，与 DRAFT
  `pair-loop-context-split` 无冲突）。`ctx.xcart` 仍是**唯一**的实验室系坐标存储（02 已把连续驱动
  的 abort 路径补齐）；`ctx.fixedatom` 即锚定位集，不复制（`pack_context.rs:604` 的
  `debug_assert_atom_props_sync` 就是这份镜像曾经漂移过的证据）；`ctx.comptype` 也只读不抄。
- **没有** `topology` 字段（评审 🔴 的落实）：本链内既没有它的生产者也没有它的消费者——01 明确把
  系统级 `Topology::tile` 推迟给 `dg-refine`，05 的"② `build_context` + `PackState` 组装"无从填它，
  04/05/06 也从不读它。加进来只会把 01 已推迟的工作拖回来（law § 3）。等 `dg-refine` 同时带来
  生产者（`Topology::tile`）与消费者时再加字段，那是一次加法扩展。
- **没有** per-atom `placed: Vec<bool>`：`OverlapField::is_placed`（`src/grow/field.rs:109`，
  由 insert/retract 增量维护、`n_placed` 计数）已经是这个 crate 的**活**放置位集，再造一份就是同一
  谓词在一次运行内的第二个真相（law § 9）。链式检查只需要状态级的形状标记
  `Placed { None | All }`，由阶段的 `guarantees()` 推进（04 / 05）。
- `Placed` 定义在这里而不是 `stage.rs`：它描述的是状态的形状，`stage.rs`（04）从 `context` import
  它——方向与 `stage.rs` 只 import `context` / `handler` / `target` / `error` 一致。
- `rigid` 是**非可选**槽（不是 `Option`）：`PackState::new(ctx, nmol)` 直接建
  `RigidView::fresh(nmol)`（`placed` 同时置 `Placed::None`）。这样"槽里到底有没有视图"这个非法状态
  不可表示（law § 8），也正是 05 得以无分支地从槽里取 `Placements` 的前提。
- 访问器：`ctx()` / `ctx_mut()` / `placed()` / `set_placed()` / `rigid()` /
  `rigid_split_mut() -> (&mut PackContext, &mut RigidView)` / `into_parts() -> (PackContext, RigidView)`。
- `invalidate_geometry_cache(&mut self)`：转调 `PackContext::invalidate_geometry_cache`
  （`pack_context.rs:529`，已存在，不重实现）。放在这里的另一个理由是 `pack_context.rs`（1042 行）
  本链一行不加。**它的调用方是 05 的阶段边界**——管线在每个阶段 `run` 之前调
  `state.invalidate_geometry_cache()`，使多阶段拼写与 `seeded_from` 拼写都从冷缓存出发（见 05 的
  几何缓存一节），因此它不是无人调用的转发器（law § 3）。
- `evaluate_unscaled(&mut self, x: &[F]) -> (F, F, F)`：**本链唯一的未缩放裁决原语**。

**未缩放裁决的合并**（本 spec 的实质工作）。今天有两套写法：

1. `src/gencan/phases.rs:33 evaluate_unscaled`——把 `radius_ini` 换进 `sys.radius`，`FOnly` 评估，
   再从 `sys.work.radiuswork` 换回；
2. `src/grow/driver.rs:381-383` 与 `src/grow/lattice/mod.rs:295-297`——`sys.scale = 1.0;
   sys.scale2 = 0.01;` 后 `FOnly` 评估。

等价性论证（逐条对源核实）：

- `sys.radius` 只在 GENCAN 路径被缩放（`src/gencan/phases.rs:189,283`、`src/movebad.rs:42`，后者在
  `:168` 还原），生长路径上 `radius == radius_ini`，半径互换是空操作。`work.radiuswork` 始终按
  `WorkBuffers::new(ntotat)`（`pack_context.rs:375`）足额分配，互换不会越界。
- `sys.scale` / `sys.scale2` 的赋值点**共三处**，不是两处：`src/grow/driver.rs:381-382`、
  `src/grow/lattice/mod.rs:295-296`，以及 `src/initial.rs:329-330`（Packmol `initial.f90:50-51`
  的移植）。三处写的都是构造器默认值（`1.0` / `numerics::DEFAULT_SCALE2 == 0.01`，
  `pack_context.rs:370-371` 的结构体字面量）。因此在任何调用点上，这两行都是把默认值写回去的
  空操作。（`src/objective.rs:414,622` 的 `let scale2 = sys.scale2;` 是读，不是写。）
- 半径互换只移动 `f_total`，不移动 `fdist` / `frest`：`fdist` 在 `src/objective.rs:292-300`
  **无条件**由 `radius_ini` 计算（在 `if overlap` 块之外），`frest` 来自约束项，它们读
  `scale` / `scale2` 但不读半径。这一条写进模块文档，让
  `pack_state_regression_unscaled_verdict_golden` 断言的是"对的理由下的那个三元组"。

合并后的单一实现按下列顺序，并且**对称地存—还两组字段**：

    存 scale / scale2 → 置 scale = 1.0, scale2 = DEFAULT_SCALE2
      → 存 radius 到 work.radiuswork，换入 radius_ini
        → EvalMode::FOnly
      → 从 work.radiuswork 还原 radius
    → 还原 scale / scale2

对称性是必要的：`scale` / `scale2` 是被公开重导出的 `PackContext` 的 `pub` 字段，而本函数的
rustdoc 承诺"返回时恢复调用方的取值"；只还 `radius` 不还 `scale` 会让承诺对一半字段撒谎
（law § 5 / § 8）。今天三处调用点的还原值都等于写入值，所以对称化本身也是逐位空操作。

实现搬到 `src/context/pack_state.rs`，以
`pub(crate) fn evaluate_unscaled(ctx: &mut PackContext, x: &[F]) -> (F, F, F)`
自由函数形式提供给持 `&mut PackContext` 的既有调用点（`gencan/phases.rs` 热路径不变形），
`PackState::evaluate_unscaled` 是它的薄转调；一个实现体，两种拼写，不是两份权威。
`gencan/phases.rs:33` 的定义删除——`pipeline/`（05）不得 import `gencan/`，这是它必须搬家的原因。

**这是一次公开面移除，在此明账（law § 10）**：`evaluate_unscaled` 今天是 `pub mod gencan`
（`src/lib.rs:99`）下 `pub mod phases` 里的 `pub fn`（`src/gencan/phases.rs:33`），所以
`molpack::gencan::phases::evaluate_unscaled` 是一条已发布路径；本 spec 把它重新安家为
`src/context/pack_state.rs` 的 `pub(crate)`，即把这条路径从公开面上撤下。`docs/architecture.md:37`
与 `docs/extending.md:460` 今天仍把它写在 `phases.rs`，两处由 07 的文档任务跟随更新（已列入 07 的
Files 与 docs 任务）。

**earn-complexity 的明账**（law § 3）：`PackState` 本身在本 spec 内除测试外没有调用方，其直接调用方
是同一条已批准链里紧随其后的 `stage-pipeline-04-stage`（seam 换成 `&mut PackState`）与
`-05-pipeline`（链式检查读 `placed`、出口读 `rigid`、阶段边界调 `invalidate_geometry_cache`）。本
spec 独立交付的价值是未缩放裁决原语的去重（三处副本、两套习语 → 一处一套）。若 04 / 05 被取消，
`PackState` 应随之回退；因为它在本 spec 内是 `pub(crate)`，这条移除路径不触及任何已发布的公开符号。

### Reuse decision

- `gencan/phases.rs:33 evaluate_unscaled` — **generalize**：搬到 `context/pack_state.rs` 成为唯一实现，`gencan/` 的定义删除（公开面移除已在上文明账）。
- `grow/driver.rs:380-386` / `grow/lattice/mod.rs:296-300` 的 scale 复位习语 — **generalize**：被同一实现吸收（等价性证明见上）。
- `pack_context.rs:529 invalidate_geometry_cache` — **reuse**：`PackState` 转调，不重实现；调用方是 05 的阶段边界。
- `pack_context.rs:209 fixedatom` / `:211 comptype` — **reuse**：经 `PackState` 只读，绝不复制。
- `grow/field.rs:109 OverlapField::is_placed` — **reuse**：保持生长内唯一的活放置位集；`PackState` 只带形状标记 `Placed`，不带位集。
- `context/rigid_view.rs::RigidView`（02）— **reuse**：`PackState` 只持槽，不扩展它的公开面。
- `context/state.rs:7 RuntimeState` — **无重叠**：借用只读遥测视图 vs 拥有型运行状态；文件命名理由成立。
- `assemble.rs:271 offset_index_column` / `Topology::tile` — **new — 本 spec 不引入 `topology` 字段**：链内无生产者也无消费者，随 `dg-refine` 一起落地（law § 3）。
- `entry/mod.rs:341 fdist: sys.fdist` — 由 05 裁决（裁决口径：读状态，不额外评估）。

## Files to create or modify

- `src/context/pack_state.rs` (new)
- `src/context/mod.rs`
- `src/gencan/phases.rs`
- `src/grow/driver.rs`
- `src/grow/lattice/mod.rs`
- `tests/context_pack_state.rs` (new)

## Tasks

- [ ] Write failing tests in `tests/context_pack_state.rs`：`PackState` 的包裹语义（`ctx()` 读到的 `fixedatom` / `comptype` 与源 `PackContext` 同一份，无副本）、`placed` 的读写、新建状态上 `rigid().nmol()` 等于构造参数且 `placed() == Placed::None`、`rigid_split_mut` 的双可变借用、`invalidate_geometry_cache` 转调生效、`evaluate_unscaled` 在缩放与未缩放半径下的返回与**存—还对称性**；确认 RED
- [ ] Add `src/context/pack_state.rs`：`Placed`、`PackState`（均 `pub(crate)`）与访问器 `ctx` / `ctx_mut` / `placed` / `set_placed` / `rigid` / `rigid_split_mut` / `into_parts` / `invalidate_geometry_cache`；在 `src/context/mod.rs` 声明模块（`src/lib.rs` 不动）
- [ ] Move the unscaled-verdict primitive into `src/context/pack_state.rs` as the single home（自由函数 + `PackState` 薄转调），对称存还 `scale` / `scale2` 与 `radius`，删除 `src/gencan/phases.rs:33` 的定义
- [ ] Switch the three call sites：`src/gencan/phases.rs` 的热路径调用、`src/grow/driver.rs:380-386`、`src/grow/lattice/mod.rs:296-300`（两处 scale 复位行被吸收）
- [ ] Add module docs to `src/context/pack_state.rs`：三个 `scale` / `scale2` 写点（含 `src/initial.rs:329-330`）、对称存还的理由、"`fdist` 由 `radius_ini` 无条件计算（`src/objective.rs:292-300`）、`frest` 读 `scale` / `scale2` 但不读半径，故半径互换只移动 `f_total`"，以及 `molpack::gencan::phases::evaluate_unscaled` 这条公开路径被撤下的明账
- [ ] Add regression scenario `pack_state_regression_unscaled_verdict_golden` to `tests/context_pack_state.rs`（固定 6 二聚体 / 20 Å 盒 / seed 7 夹具，硬编码 `(f_total, fdist, frest)` 金标，容差 1e-12）
- [ ] Run full check + test suite

## Testing strategy

- 归属：`PackState` 的契约归 `tests/context_pack_state.rs`（law § 11；夹具用
  `PackContext::new(ntotat, nmol, ntype)` 自造，与 `tests/geometry_cache.rs:28` 同法）。单测门：
  `cargo test -p molcrafts-molpack --lib --tests -- pack_state`。
- Happy path：包裹后各访问器可达；`Placed` 状态迁移；`rigid_split_mut` 的双可变借用编译且互不别名。
- Edge cases：`ntotat == 0` 的退化上下文；`evaluate_unscaled` 在 `radius != radius_ini` 时正确还原
  `sys.radius`（互换不泄漏）；在 `scale != 1.0` / `scale2 != DEFAULT_SCALE2` 的人为状态上调用后，
  两个字段被**还原成调用前的值**（对称性直测，今天的实现会红）。
- **值等价（本 spec 的核心）**：同一夹具在改动前后，`evaluate_unscaled` 返回的
  `(f_total, fdist, frest)` 与调用后 `sys.radius` 的内容逐位相同（`to_bits()`）——GENCAN 夹具
  （半径已被 discale 缩放）与生长夹具（`radius == radius_ini`）各一。
- Regression scenario：`pack_state_regression_unscaled_verdict_golden`，硬编码字面量金标；
  对应 `type: runtime` 验收项。
- 既有逐位门：`tests/gencan.rs`、`tests/grow.rs` 确定性组、`tests/packer.rs`、
  `examples_batch`（release，`--ignored`，五例）。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- `topology` 字段：本链无生产者也无消费者，随 `dg-refine` 的 `Topology::tile` 一起落地（law § 3）。
- `scale`（半径阶梯）与 `bond_tolerance` 字段：同上。
- per-atom `placed` 位集——永久不做，`OverlapField::is_placed` 是它的家。
- `PackState` / `Placed` 升为 `pub`（04）。
- `docs/architecture.md:37` / `docs/extending.md:460` 里 `evaluate_unscaled` 位置的更新（07 的文档任务）。
- 把 `PackContext` 的字段物理搬进 `PackState`；若将来要做，必须与 `pair-loop-context-split` 合成
  一份 objective / context 重组 spec。
- seam 签名改动（04）与生命周期体搬迁（05）。
