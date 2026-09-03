---
title: stage-pipeline-02-view — RigidView：刚体自由度的唯一归属
status: approved
created: 2026-09-02
chain: stage-pipeline（02 of 7）
---

# stage-pipeline-02-view

## Summary

把"刚体放置向量"这件事收进一个类型 `RigidView`（`src/context/rigid_view.rs`）：它吸收今天的
`PlacementsMut` 访问器、`initial::init_xcart_from_x` 的两个方向、`GenCanPack::prepare` 的种子注入、
以及两个生长驱动里逐字重复的写回块。前置一步把 `xcart` 变成实验室系坐标的唯一家——连续驱动的 abort
补齐路径今天不回写 `xcart`，本 spec 补上，与 lattice 驱动（`src/grow/lattice/mod.rs:262`）取齐。
`RigidView` 只装刚体自由度本身，**不带**任何"这份放置是被喂进来的"标记——那件事的家是 03 的
`PackState.placed`；`GencanSettings.push_off` 本 spec **一行不动**（仍由 `src/gencan/entry.rs:213`
从 `seed_placements.is_some()` 推出），05 用 `Placed` 取代它。生命周期第 ⑤ 段因此不再需要
`entry/mod.rs` import `crate::initial`。三个入口的坐标逐位不变。

## Design

**层次**：共享状态底座（`src/context/`）。其余被触及的文件都是被搬迁符号自己的调用点。

**为什么不是 `src/gencan/view.rs`（父 spec 的位置）**：`PackState`（03）要持有 `rigid: RigidView`，
那会造成 `crate::context → crate::gencan` 的 import 边。`context/rigid_view.rs` 是同时满足
law § 4 / § 6 的唯一位置；`init_xcart_from_x` 的实现只用到 `crate::euler`（叶子）与 `PackContext`
字段，搬到这里没有任何新依赖。

**`context` 不得认识 `entry`（评审 🚨 的直接落实）**：`Placements`（`src/entry/result.rs:18`）是
`entry` 层的出口快照类型；`src/context/` 今天对 `crate::entry` 零引用，而 `src/entry/mod.rs:23`
已经 import `crate::context`。因此 `install_seed` **不收 `Placements`**，只收两片裸数据：

    pub fn install_seed(x: &[F], coor: &[[F; 3]], ctx: &mut PackContext) -> Self

拆包在调用方 `src/gencan/entry.rs` 完成（`RigidView::install_seed(seed.rigid.as_slice(),
&seed.coor, sys)`）。方向保持 `entry → context` 单向，acceptance 用 grep 钉住。

**先把 `xcart` 变成唯一的家（评审 🔴 的直接落实）**：连续驱动今天只在**回合末**把
`chain.coords` 同步进 `sys.xcart`（`src/grow/driver.rs:318-321`）；abort 补齐循环
（`:352-357`）用 `force_place` 改了 `chain.coords` 却不再同步，所以今天的写回块（`:360-378`）
必须从 `chain.coords` 取值才正确。lattice 驱动没有这个洞（`src/grow/lattice/mod.rs:262` 在
abort 循环内就写了 `xcart`）。本 spec 的**第一步**是在连续驱动的 abort 补齐循环内补上同一行同步，
使 `ctx.xcart` 在写回发生前就是权威（law § 9）；只有在这一步之后，`capture_from_xcart` 从
`ctx.xcart` 读取才与今天逐位相同。这一步单独可验证：先在 `tests/grow.rs` 落一条
`grow_abort_writeback_golden`（用 `EarlyStopHandler` 触发 abort，把今天构建下的
`PackResult::positions()` 硬编码为字面量），它在同步前后都必须绿。

**类型**（`src/context/rigid_view.rs`，预算 ≤ 300 行）：

    pub struct RigidView { x: Vec<F>, nmol: usize }   // Debug + Clone

- **布局契约**：`x` 是今天的扁平向量——先 COM 块后 Euler 块，每分子各 3 个（`PlacementsMut` 的
  文档契约逐字继承）。访问器 `nmol` / `com` / `set_com` / `euler` / `set_euler` / `as_slice` /
  `as_mut_slice` 原样搬入；`RigidView::fresh(nmol)` 取代 `vec![0.0; 6*nmol]` + `PlacementsMut::new`，
  长度断言随之内化（非法状态不可表示，law § 8）。
- **没有 `seeded` 旗标**：视图只知道自己持有的 6N 个自由度。"这些放置是不是必须被后续阶段接续"
  是**状态**级的事实，家在 03 的 `PackState.placed`（law § 9）。若视图也存一份，GENCAN 阶段用
  `set_com` / `as_mut_slice` 写完解之后不会去更新它，两份表示当场分叉——这正是评审 r2 指出的
  第二个家。rustdoc 写明这条缺席及其理由，免得后人"顺手补上"。
- `write_xcart(&self, ctx: &mut PackContext)` = `src/initial.rs` 的 `init_xcart_from_x` 逐字搬迁
  （`xcart = com + R(euler)·coor`）。**两个方向都归它**：第 ⑤ 段的出站重建
  （`src/entry/mod.rs:299`）与 push-off 的入站重建（`src/gencan/solver.rs:161`）是同一个函数的
  两次调用；`src/initial.rs::init_xcart_from_x` 删除，无别名。
- `install_seed(x: &[F], coor: &[[F;3]], ctx: &mut PackContext) -> Self` = `src/gencan/entry.rs:193-195`
  的种子注入（`ctx.coor[..n].copy_from_slice(coor)`，视图的 `x` 取自 `x`）。这是**构造形**：种子
  来时还没有视图，因此它**返回**视图、由调用方写进自己的槽。本 spec 里调用方是
  `GenCanPack::prepare`，它把返回的视图写进本次 run 的 `x` / 视图缓冲（与今天
  `x.copy_from_slice(&seed.x)` 同位）；从 03 起槽由 `PackState` 持有，写法统一为
  `*view = RigidView::install_seed(..)`。网格与 simbox 的安装（`:181-192`）留在原处，本 spec 不动。
- `capture_from_xcart(&mut self, ctx: &mut PackContext)` = `src/grow/driver.rs:360-378` 的写回契约
  提为唯一实现，`src/grow/lattice/mod.rs:274-292` 的逐字副本删除、改为调用。**填充形**：填的是
  已存在的视图（03 的 `PackState` 持槽，04 把 `&mut RigidView` 交给阶段），不返回新值。`ctx` 取
  `&mut` 是因为写回契约的后半段要把居中构象存进 `ctx.coor`（正是 `write_xcart` 的逆），把两半拆成
  两个方法会造出一个调用方能忘的步骤（law § 8）。它按 `ctx` 的 `(idfirst, natoms, nmols)` 布局
  枚举拷贝，因此 `x` 的分子下标恒等于 xcart 布局下标。

**`Placements` 的形状**（`src/entry/result.rs:18`）：`x: Vec<F>` 字段换成 `rigid: RigidView`，
`coor` / `copy_atoms` / `cell` 不动。运行期 `coor` 的唯一权威仍是 `ctx.coor`；`Placements` 是出口
快照（law § 9 明确允许的表示副本，不是第二份可变真相）。方向是 `entry → context`，不新增反向边。

**seam 签名**：`Solver::solve` 的 `x: PlacementsMut<'_>` 改为 `x: &mut RigidView`，
`pub struct PlacementsMut` 删除（不留别名，law § 1）。`Solver` / `SolveOutcome` 的名字本 spec 不动。

**不做的事（对上一稿的撤回）**：上一稿声称 lattice abort 路径存在 `mol` / `base` 错配。复核
`src/grow/lattice/mod.rs:156-208`（`done` 按布局序推入，`mol` 与同一次 `(itype, imol)` 走查同步
递增）与 `:230-232`（唯一的 `break 'outer` 在 `done.push` 之后）后，abort 循环从
`m = done.len()` 重启并按同一布局序跳过已完成的 base，`mol` 与布局下标始终一致——该缺陷不可达。
本 spec 撤回该断言与它的验收项；若将来出现反例，它属于 `/mol:debug`，不折进重构（law § 10）。

### Reuse decision

- `initial.rs init_xcart_from_x`（调用点 `entry/mod.rs:299`、`gencan/solver.rs:161`）— **reuse**：函数体逐字成为 `RigidView::write_xcart`，两个方向都归它，原函数删除。
- `gencan/entry.rs:170-197 prepare` 的种子注入 — **reuse**：逐字成为 `RigidView::install_seed`（去 `Placements` 参数），`seeded_run_contract` 的逐位接续因此不变。
- `grow/driver.rs:360-378` 写回契约 — **reuse**：提为 `RigidView::capture_from_xcart`，两个生长驱动都调用它，不重新推导。
- `grow/lattice/mod.rs:262` 的 abort 内 xcart 回写 — **pattern**：连续驱动照抄这一行，使 `xcart` 成为唯一的家。
- `gencan/solver.rs:48-51 GencanSettings.push_off` — **reuse — 本 spec 一行不动**：它仍由 `src/gencan/entry.rs:213` 从 `seed_placements.is_some()` 推出；05 用状态词汇 `PackState.placed == Placed::All` 取代它并删除该字段。视图不承接这份事实。
- `src/solver.rs:63 PlacementsMut` — **generalize**：其访问器成为 `RigidView` 的公开面，类型名删除。
- `assemble.rs` 的平铺函数、`OverlapField::is_placed`、`invalidate_geometry_cache`、`evaluate_unscaled` — 与本 spec 无接触，分别由 01 / 03 裁决。

## Files to create or modify

- `src/context/rigid_view.rs` (new)
- `src/context/mod.rs`
- `src/solver.rs`
- `src/initial.rs`
- `src/entry/mod.rs`
- `src/entry/result.rs`
- `src/gencan/entry.rs`
- `src/gencan/solver.rs`
- `src/grow/driver.rs`
- `src/grow/lattice/mod.rs`
- `src/lib.rs`
- `tests/context_rigid_view.rs` (new)
- `tests/grow.rs`

## Tasks

- [ ] Write failing tests in `tests/context_rigid_view.rs`：布局契约（从 `tests/grow.rs:197,231` 的 `placements_view_layout` / `placements_view_rejects_bad_len` 迁入并改名）、`write_xcart` 对已知 `(com, euler, coor)` 的结果、`install_seed` 的逐位注入（返回视图的 `x` 与 `ctx.coor` 均等于输入）、`capture_from_xcart` 的 COM/居中/`euler = 0`、`fresh(nmol)` 的长度与 `nmol()`；确认 RED
- [ ] Pin the continuum abort writeback with `grow_abort_writeback_golden` in `tests/grow.rs`（`EarlyStopHandler` 触发 abort，硬编码今天构建下的 `PackResult::positions()` 字面量），随后在 `src/grow/driver.rs:352-357` 的 abort 补齐循环内同步 `sys.xcart = chain.coords`（照 `src/grow/lattice/mod.rs:262`），该金标改动前后都必须绿
- [ ] Add `src/context/rigid_view.rs`：`RigidView { x, nmol }` + `fresh` / 访问器 / `write_xcart` / `install_seed` / `capture_from_xcart`，rustdoc 写明视图**不**持任何 seeded 标记及其理由（该事实归 03 的 `PackState.placed`）；在 `src/context/mod.rs` 与 `src/lib.rs` 导出
- [ ] Replace `PlacementsMut` with `&mut RigidView` on `Solver::solve` in `src/solver.rs` 并更新三个实现者的签名（`src/gencan/solver.rs`、`src/grow/driver.rs`、`src/grow/lattice/mod.rs`）
- [ ] Move `init_xcart_from_x` out of `src/initial.rs`：删除原函数，`src/entry/mod.rs` 第 ⑤ 段与 `src/gencan/solver.rs:161` 改调 `RigidView::write_xcart`，`src/entry/mod.rs` 移除 `use crate::initial`
- [ ] Replace both growth writeback blocks with `RigidView::capture_from_xcart`（`src/grow/driver.rs:360-378`、`src/grow/lattice/mod.rs:274-292`），改为从 `ctx.xcart` 读取
- [ ] Switch the seed injection in `src/gencan/entry.rs::prepare` to `RigidView::install_seed`（在该处拆包 `Placements`，把返回的视图写进本次 run 的槽），并让 `src/entry/result.rs::Placements` 改持 `rigid: RigidView`（`GencanSettings.push_off` 本 spec 不动）
- [ ] Add regression scenario `rigid_view_regression_xcart_and_capture_golden` to `tests/context_rigid_view.rs`（硬编码 `(com, euler, coor)` → `xcart` 与反向 capture 的金标，容差 1e-12）
- [ ] Verify the bitwise gates：`tests/packer.rs::seeded_run_contract`、`free_chain_push_off_deterministic`、`tests/grow.rs` 确定性组，以及 `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored` 五例
- [ ] Run full check + test suite

## Testing strategy

- 归属：`RigidView` 的契约归 `tests/context_rigid_view.rs`（law § 11，一个文件对一个模块的公开
  API；夹具用 `PackContext::new(ntotat, nmol, ntype)` 自造，与 `tests/geometry_cache.rs:28` 同法）。
  单测门：`cargo test -p molcrafts-molpack --lib --tests -- rigid_view`。
- Happy path：布局读写、`write_xcart` 的旋转合成、`install_seed` 后视图的 `x` 与 `ctx.coor` 逐位等于
  输入、`capture_from_xcart` 的三条写回不变量（COM = 质心 / `coor` 居中 / `euler = 0`）。
- Edge cases：`fresh(0)`；长度不匹配 panic；多类型多拷贝布局下 `x` 的分子下标与
  `(idfirst, natoms, nmols)` 一致；单原子拷贝的 COM = 该原子。
- xcart 单一归属：`tests/grow.rs::grow_abort_writeback_golden` 是**特征化金标**——它先于同步落地并
  在同步后仍绿，证明 abort 路径改读 `ctx.xcart` 没有移动任何数值；配合既有的
  `grow_abort_keeps_bonded_geometry`（`tests/grow.rs:1398`）不动断言。
- Regression scenario：`rigid_view_regression_xcart_and_capture_golden`，浮点金标以 1e-12 断言，
  注释记录捕获自本 spec 前的构建；对应 `type: runtime` 验收项。
- 跨模块逐位门（既有回归，不是本模块单测）：`tests/packer.rs::seeded_run_contract`、
  `free_chain_push_off_deterministic`、`tests/grow.rs` 三个确定性测试、`examples_batch` 五例。
- 既有 RED `tests/grow.rs::grow_cg_kremer_grest_c_inf` 不在门内，不得削弱或跳过。
- 无物理新增，故无 Domain basis 一节。

## Out of scope

- `Solver` → `Stage` 改名与 `requires` / `guarantees` 的引入（04）。
- `PackState`（03）；本 spec 只让 `RigidView` 独立成立。
- `GencanSettings.push_off` 的删除与 push-off 判据的状态化（05）。
- `prepare()` 的删除与网格安装的搬迁（05）。
- `initial ↔ gencan` 的既有环（`initial.rs:24` / `gencan/solver.rs:18-19`）——记录不动。
- lattice abort 路径的策略与分子下标：复核后不存在错配，本 spec 不动（见 Design 末段）。
- 连续驱动 abort 补齐的**策略**（是否补齐、如何补齐）不变，只补 `xcart` 同步。
