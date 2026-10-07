# .claude/notes/notes.md — molpack evolving decisions

Use `/mol:note <decision>` to add or supersede an entry. Rules that become absolute move to `law.md` (one line each in CLAUDE.md); this file keeps the reasoning and the history.
Format per entry:

```
## YYYY-MM-DD — short title
<rule>
**Why:** <reason — incident, deadline, or stakeholder ask>
**How to apply:** <when / where this kicks in>
```

## 2026-10-07 — molrs S1–S4 命名波（`ir/p-registry` @ 166cd540）与容器词清扫

- **molrs 路径**：`molrs::{store, system, spatial, units}` → `molrs::core`（Python `molrs.core`）；`op` 扁平（`op::rigid::nerf` → `op::place_from_internal_coords`，`about`/`apply` → `rotation_about`/`transform_point`）；`LBFGS` → `Lbfgs` + `LbfgsSettings`；`Optimizer::run` → `minimize`；`OptReport` → `OptimizationReport`（多 `final_grad_rms`）；密度盒用 `core::constants::ANGSTROM3_PER_CM3`。
- **取代 K1**：molrs 删了按扩展名分派的 `io::read_frame` / `write_frame`，`.inp` 的 `filetype` 词表回到 molpack：`script::StructureFormat`（原 12 种格式，名或扩展名；`resolve` 恒编译，`read`/`write` 走 molrs 各格式自己的门，`io` 门控）。CLI、`Script::build` 与 wheel 默认加载器（同一 `resolve` + 对应 `molrs.io.read_<fmt>`）共用这一张表。
- **容器词与缩写**：`Handler` → `Callback`（`with_callback`、`ProgressCallback`、`LammpsLogCallback`、`EarlyStopCallback`、`XyzTrajectoryCallback`，Python 协议 `molpack.Callback`）；`entry/` 拆为 `settings.rs` / `pack_space.rs` / `state.rs`；`*/entry.rs` → `gencan_pack.rs` / `cbmc_grow.rs` / `lattice_grow.rs`；`testutil` → `test_fixtures`；绑定 `entry.rs` → `packing_methods.rs`，`handler.rs` → `callback.rs`，`result.rs` → `state.rs`，`types.rs` 并入 `target.rs`；`GenCanPack` → `GencanPack`，`InvalidPBCBox(Error)` → `InvalidPbcBox(Error)`。
**Why:** molrs 每个符号一条路径、缩写按词大写；molpack 自己的名字同一规则。
**How to apply:** 逐位对拍（对 `parity-base`）：四个 CLI 结构文件、`pack_peo grow` 与其余 stdout 逐字节一致；`adsorption.xyz` / `translocation.xyz` 与上一轮（molrs 0.16 写出器）逐字节一致，对基线的差别是写出器（extended XYZ，全精度），原子与元素相同、坐标差 ≤ 5e-5（基线 `%.4f` 的舍入）；`pack_peo lattice` stdout 只差 `GenCanPack` → `GencanPack` 一词。

## 2026-10-07 — 0.4.0：模块单一职责重构（wave K），对 molrs `ir/p-registry` @ 64afcf90

每个模块一个职责；molpack 不重复 molrs；每个公开符号只有一条路径；删别名、重复、死代码与兼容垫片（审计 `molnex/.claude/specs/module-responsibility-audit-2026-10-06.md` §4 K1–K11 与裁决 10）。
- **molrs 0.16 单一路径**：`molrs::core::{Frame, Block, Column, BlockDtype}`、`molrs::core::{Topology, BondDistanceWeights, Element, Atomistic, Atom, TopologyError}`、`molrs::core::{SimBox, Mic, BoxKind, TriMesh}`、`molrs::op::{F, Idx, …}`；`SoftLbfgs` → `LBFGS::new(Arc::new(SoftSpec::from_frame(..).potential(None)), …)`（`SoftSpec` 在 `molrs::ff::potential::soft`）。
- **K1** 删 `src/script/io.rs`（`molpack::script::{read_frame, write_frame}`）：`Script::build` 与 CLI 走 `molrs::io::{read_frame, write_frame}`；wheel 的默认加载器就是 `molrs.io.read_frame`。格式表只在 molrs（CLI 因此多读写 MOL2/GRO/CIF/POSCAR/XSF/cube/inpcrd/LAMMPS data 写出）。
- **K2** `XYZHandler` 改为 `io` 门控，快照拼成 Frame（`element`/`x,y,z`/`mol_id`，meta `step`）交给 `molrs::io::data::xyz::write_xyz_frame`；文件格式变为 molrs 的 extended XYZ（`step=N Properties=species:S:1:pos:R:3:mol_id:I:1`）。
- **K3** `assemble::topology_frame` 改为每个目标 `Frame::replicate` 后 `Frame::concat`，删掉手写的端点偏移与 `row_bases`。**K4** `molrs::core::constants::AVOGADRO`。**K5** 平面重复判断问 `HalfSpace::new(n, 0)?.repeats_along(s)`。
- **K6** grow 内坐标用 `op::vec3::{angle, dihedral}` + `op::rigid::nerf`（逐位等同旧实现）；集体约束法向用 `op::vec3::normalize`（Python 绑定同一判据）；`TorsionMcOptimizer` 的 `recenter_free` 用 `op::superpose::centroid`，退化键判据用 `normalize`。**裁决 10 保留**：`.inp` 几何约束核（`restraint/geometric`，逐位复现 Packmol `comprest`/`gwalls`）、金刚石晶格行走（`grow/lattice/saw`）、独立 RNG 流（`random.rs`）、模板中心的求和顺序（`target::geometric_center`）——各自文档写明原因。
- **K7** 删根上的 molrs 再导出（`F`、`Element`、`BondDistanceWeights`、`Optimizer`、`OptimizationReport`）与 `prelude`。**K8** 只剩三个公开命名空间：`context`（布局常量、`AtomProps`、`WorkBuffers`、`GeometryKey`）、`grow`（`GrowConfig`、`GrowError`、`LatticeConfig`、`TorsionPrior`、`AnglePrior`）、`script`；其余模块私有，根上一条路径（新增根导出 `StageProgress`、`EngineSetup`、`Restraint`、`GroupEvaluation` 与六个集体约束）。原模块文档中面向用户的内容移到 `Pipeline` / `Stage` / `Invariant` / `AtomRestraint` 的类型文档。
- **K9** `compute_f/fg/g` 改 `pub(crate)`；删 `PackContext` 的固有 `evaluate`（唯一入口 `Objective::evaluate`）、`EvalMode::RestMol`、`grow::moves::uniform`。顺手删死代码：`GencanParams::{iprint, ncomp}`、`GencanResult::{fcnt, gcnt, cgcnt}`、`CgResult::{q, iter}`、`OverlapField.n_placed` 字段、`InternalTree.n_atoms` 字段（测试访问器改为推导并 `cfg(test)`）、`ATOM_PROPS_SIZE`（改 `const _: () = assert!(..)`）。
- **K10** 唯一版本闸门是扩展导入时的 `molrs_capsule::check_abi`：删 `molpack/version.py`（`MOLRS_MINOR`、`check_molrs_version`）及其测试；`molpack.version` 由 `__init__.py` 读 wheel 元数据。唯一 CLI 是 Rust `molpack` 二进制（README 与 `docs/cli/` 呈现的就是它：`molpack mixture.inp` / stdin）；删 Python `molpack.cli`（typer 的 `molpack pack|info|version`，语法不同）、`[project.scripts]` 与 `typer` 依赖。
- **K11** 解散 `numerics.rs`：`objective_small_floor` ≡ `residual_small_floor` 合为 `pack/gencan` 的 `small_floor`，其余 floor 也归 GENCAN，`numeric_controls` 归 `search.rs`，`DEFAULT_SCALE2` 归 `context/pack_context.rs`；`grow_error.rs` → `grow/error.rs`；Python `helpers.rs` → `errors.rs`（删 `NpF = f64`，用 `molrs::op::F`），`constraint.rs` → `restraint.rs`（`tests/test_constraint.py` → `test_restraint.py`）；删 `ff` 透传特性——唯一用户是一个测试，改为 molrs dev-dependency 带 `ff`；门禁命令 `--features cli,ff` → `--features cli`。
**Why:** 所有者规则（2026-10-06）：每模块一职，molpack 用 molrs 的 API，不复制。
**How to apply:** 升 molrs 小版本的清单去掉 `version.py::MOLRS_MINOR`（文件已删；取代 0.3.0 条目中的这一项）。逐位对拍：对同一 molrs 导出，path-only 基线（a383ff6 + 路径改写）与重构版在 CLI `mixture`/`bilayer`/`interface`/`solvprotein`、`pack_peo grow`/`lattice`、`pack_adsorption`、`pack_translocation` 上输出逐字节一致；lib 381（少的是随 `script/io.rs` 删除的读 dump 测试）+ doc 22（少的是 `prelude` 的）、clippy `--all-features -D warnings`、`--no-default-features`/`rayon`/`io` 检查、rustdoc 零警告、fmt 全绿；wheel 测试 144 通过（少的 4 个是 `test_version.py`）。

## 2026-10-06 — 0.4.0 发布线：molrs 0.16（未发布）

molpack 0.4.0 跟随 molrs **0.16** 小版本线（molpack 小版本随 molrs 小版本各升一级：0.3 ↔ 0.15，0.4 ↔ 0.16）。按 0.3.0 条目的清单同步改了 `Cargo.toml`（`version = "0.16"`）、`python/Cargo.toml`（molrs + molrs-ffi `0.16`，molpack `0.4.0`）、`pyproject.toml`（`molcrafts-molrs>=0.16.0,<0.17`；`[molpy]` extra 与 typecheck 组 `molcrafts-molpy>=0.16.0,<0.17`；tox 的两处 `(0,16)` 断言）、`version.py::MOLRS_MINOR = (0, 16)`、两个发布工作流的 `MOLRS_GIT_REF` / `mol_git_ref` = `v0.16.0`。路径仍是仓库相对的 `../molrs/molrs`（CI 布局）。
molrs 0.16 的破坏性变更（力场 IR 协议：`CompileError`、`register_kernel*` 返回 `Result`、typifier `Match.links`、`Style::category() -> &str`、谐振 K 不再减半、角度改度；`molrs.md` 不再导出 `LJCut`/`Potential`/`Potentials`）都不触及 molpack：molpack 不构造力场、不编译势函数，`ff` 只是透传 `molrs/ff`（文档只按名字提到 `LBFGS`，0.16 仍在 `molrs::optimize`）；Rust 与 Python 两侧都无需改代码。已在 CI 布局下对 molrs `ir/p-registry` @ 2bb03675（版本改为 0.16.0 的导出副本）验证：lib 382 + doc 23 测试、clippy `--all-features`、`--no-default-features` / `rayon` 检查、fmt 全绿；wheel 测试 148 通过（molrs 0.16.0 wheel + molpack 0.4.0 wheel，Python 3.12）。
Python 端 relaxer：维持 0.3.0 条目的裁决——in-loop 优化器的 Python 绑定仍是待设计项，`feat/subset-optimizer` 的 `molpack.relaxer` / `SoftOptimizer` 不移植。
**Why:** molrs 0.16 改了力场 IR 与 ABI 线（capsule 名带 `/0.16`），两个 wheel 必须同一小版本线。
**How to apply:** `v0.16.0` 标签在 molrs 远端还不存在（molrs 0.16 未发布、对应提交也未推送）；发布 molpack 0.4.0 之前先确认 molrs 打出 `v0.16.0`，否则发布工作流检出失败。本地验证用 CI 布局：sibling `molrs/` 是 0.16.0 版本的 molrs 源码树。

## 2026-10-06 — 0.3.0 发布线：molrs 0.15.0 标签；feat/subset-optimizer 由优化器接缝取代

molpack 0.3.0 对 molrs **v0.15.0 发布标签**构建（不是 molrs dev）。发布工作流（`publish-crate.yml` / `publish-pypi.yml`）检出 `v0.15.0`；`ci.yml` 在 `workflow_call` 上接 `mol_git_ref` 输入，发布时传 `v0.15.0`。molpy 不是运行时依赖，只给 `pack_peo_*.py` 用：`[molpy]` extra 钉 `molcrafts-molpy>=0.15.0,<0.16`。
`feat/subset-optimizer`（`SoftOptimizer`、`Target::with_optimizer`、`Molpack::with_subset_optimizer`）不合入：它挂在已删除的 `Molpack` 上。Rust 侧由 `GenCanPack::with_optimizer` + `OptimizeSelect::{per_copy, joint}` + molrs `SoftLbfgs`/`SoftSpec`（软 overlap + 1-2/1-3 的唯一实现在 molrs）取代；约束项改由非损害门（整体目标含约束，变差即回滚）把关。joint 契约的单测补在 `src/optimizer/mod.rs`。
**Why:** dev 原本对 molrs 发布前的 API 编译（`relation_endpoints` 两参、`Block::get_bool`、Python `Block.view`），对 v0.15.0 标签编译失败；tag 触发的 CI 会把 molpack 标签名当作 molrs/molpy 的 ref。
**How to apply:** 升 molrs 小版本时同时改 `Cargo.toml` / `python/Cargo.toml` / `pyproject.toml` / `version.py::MOLRS_MINOR` / 两个发布工作流的标签（含 `jobs.ci.with.mol_git_ref`）。Python 端的 in-loop 优化器绑定（0.2.0 的 `TorsionMcRelaxer` / `LBFGSRelaxer` / `molpack.relaxer`）在 0.3.0 中没有替代品，属于待设计项，不要用旧 `Molpack` 形状补回。

## 2026-10-02 — 不设 WalkGrow；晶格驱动就是自回避随机行走

不增加理想链入口 `WalkGrow`。自回避随机行走是 `LatticeGrow` 的职责。CBMC 与晶格继续是两个驱动，不合成 `Grow<S, X>`。
**Why:** 操作者 2026-10-02：理想链入口没有单独的职责，晶格驱动本身就是自回避随机行走。
**How to apply:** 熔体生成走晶格。不要为无排除体积的连续行走再开一个入口。`dg-refine` 的重叠夹具用人工重叠放置，不依赖 `WalkGrow`。
**Supersedes:** grow-axes 里的 `WalkGrow` 预设，以及 2026-09-02「熔体主路径改为理想链生成」中的理想链入口。

## 2026-09-29 — 清理：删死代码、优化器接缝去门控、GenCanPack 默认早停

- **删除**：`src/cases.rs`（`ExampleCase`、`build_targets`、`example_dir_from_manifest`、`render_inp_script`）；校验模块 `validation`（`validate_from_targets`、`ValidationReport`、`ViolationMetrics`、`MOLPACK_DEBUG_VALIDATION`）——裁决只读 `State::fdist`/`frest`，阶段守卫是 `RestraintsSatisfied`；`src/system/state.rs`（`RuntimeState`/`RuntimeStateMut`、`runtime()`/`runtime_mut()`）与 `src/system/model.rs`（`ModelData`、`model()`）——调用方直接读 `PackContext`/`PackState`；`template::frame_positions` 与 `FramePositionsError`——改为 `molrs::core::Frame::coords` + `template::coord_rows`；`NullHandler`、`Handler::on_inner_iter`、`StepReport.relaxer_acceptance`；`PackContext` 的 `RestraintRef` 别名（就是 `usize`）。单测共享夹具在 `src/testutil.rs`（`cfg(test)`）。
- **优化器接缝去门控**：`src/optimizer/`、`GenCanPack::with_optimizer`、`OptimizeSelect`、`TorsionMcOptimizer` 始终编译；`ff = ["molrs/ff"]` 只是透传，molpack 自身不因它多编译任何东西。「relaxer」概念整体退场，文档页为 `docs/rust/handlers-optimizers.md`。
- **默认早停**：`GenCanPack` 默认装 `EarlyStopHandler`（`src/handler.rs`，与 Packmol 对齐：radscale = 1 时 10 轮内 bestf 改善 < 10 % 即结束该阶段）；`with_early_stop(None)` 关闭、`with_early_stop(handler)` 替换；`GenCanPack::default_max_loops(ntype) = 200 * ntype`。

**Why:** 删除项在 `tests/`/`benches/`/`regressions/` 退场（2026-09-20）后已无消费者；接缝只依赖 molrs 核心的 `Optimizer`，`ff` 门控没有理由；早停对齐 Packmol 的平台期收手，而不是跑满 `max_loops`。
**How to apply:** 新代码不要再引用上述名字；判收敛读 `State`，读模板坐标用 `Frame::coords` + `coord_rows`；想跑满预算显式 `with_early_stop(None)`。

## 2026-10-01 — pack 一族、grow 一族；D-07 清

刚性放置（`initial` / `movebad` / `gencan`）在 `src/pack/`，对箱外 `pub(crate)`。生长仍是 `src/grow/`，晶格行走是它的同级子模块，不是 CBMC 驱动的孩子。`euler` 与可选的循环内优化器留在族外。网格安装在 `system::grid`（`CellGrid::for_cutoff_capped`），生长不再依赖 GENCAN 驱动。`GrowError` 在箱根叶子 `grow_error`（晶格变体留在同一个枚举上）。`EvalMode` / `EvalOutput` 在 `eval`；`PackContext::evaluate` 的实现在 `objective`；`Constraints` 删除。`StageOutcome` 在 `outcome`。`Target::fixed_from` 收 `&Frame`。阶段类型名是 `GenCanStage`（`NAME` 仍是 `"gencan"`）。箱根不再重导出 `LatticeConfig`。

**D-07 清**：(1) 单点卷回已是 `SimBox::wrap_row`；(2) 网格封顶是 `CellGrid::for_cutoff_capped`；(3) `topology_for_growth` 映射 `TopologyError`（`MissingBondEndpoint` / `BondOutOfRange { row, atom, n }`），不再自扫键。`Block` 读列走 `Column::as_*`。
**Why:** 审查要求一族一个模块，并且和本地 molrs 0.15 对齐。环（initial↔gencan↔movebad、objective↔constraints、handler↔stage、target↔entry、error↔grow）是违规，不是风格。
**How to apply:** 新的刚性驱动放进 `pack/`；生长的拒绝只加在 `GrowError` 上；网格覆盖只读 `system::grid::coverage_radmax`（直径）。不要打开 molrs `builder`，也不要把 CBMC 和晶格并成一个泛型驱动。

## 2026-09-29 — 对齐 molrs 0.15 定稿

molpack 以 `default-features = false` 依赖 molrs / molrs-ffi（molrs 的 default 是 `full` 全家桶，此前 molpack 的 `io`/`ff`/`rayon` 特性只管自己的代码，管不住依赖）。同批清掉的 molpack 侧问题：优化器环境的最小像改用 `SimBox::mic()`（原实现只要任一轴周期就三轴全卷，slab 盒错）；`mol_id` 统一 1 起；`assemble` 改用 molrs `Block::select_rows`/`merge`/`Column::resize`，回放 schema 的全部关系块（含 `pairs`/`exclusions`），模板列 dtype 冲突为具名错误 `PackError::TemplateColumns`（运行前检查）；未分类键→可旋转单键的策略只在 `template::rotatable_bonds` 一处；`TorsionMcOptimizer` 排除表来自 molrs `Topology::exclusions`，`with_special_bonds` 可改；晶胞只剩 `SimBox`（`declared_cell -> Option<SimBox>`，删 `CellDeclaration` 与 `ResolvedSpace.pbc`，`CellDecl` 只是构建器的未校验输入）；氢的识别是 `Target::with_hydrogens`（默认元素 `H`）；删 `src/frame.rs` 与 `PackContext.frame`。

**D-07** 已于 2026-10-01 清偿，见上条。
**Why:** law「一事一家」「依赖随策略」。
**How to apply:** 见 2026-10-01 条。

## 2026-09-20 — 审查修复：一套 cell list、一处网格覆盖、退休无生产者的声明面

- **D-02 清（是真 bug，不是理论问题）**：`install_resolved_cell` 用 `max(radius_ini)`、`initial()` 用 `2·max(radius_ini)` 作网格 `radmax`。配对核的作用距离是 `radius_i + radius_j`（已乘 `discale`），即 `2·discale·R`，而 ±1 stencil 只保证找到相距 < cell_side 的对——前者的 cell_side = `1.01·discale·R` 只有需要的一半，**生长 / lattice / 接续 GENCAN 三个阶段前奏真的漏配对**（新 RED 测试 `initial::grid_coverage_tests` 证明：2.15 Å 的重叠对在 2.2 Å 接触下被判为零惩罚）。两处统一到 `initial::coverage_radmax`（直径，读 `radius_ini`）。
- **一套 cell list**：`grow/field.rs::OverlapField` 自建的 `cell_index` / `unflat` / `shift` 27-邻域与最小镜像全部删除，改用 `molrs::core::CellGrid`（周期轴 wrap、自由轴 clamp、小网格去重）与 `SimBox::mic()`。molpack 只保留**占用**（链表 + 空 cell 簿记），不再持第二套格点划分。净 −57 行。
- **热路径分配**：`objective::accumulate_collective_fg` 每次求值分配的 `spans` 挪进 `WorkBuffers::collective_spans`（take / put back，与 `collective` 同款）。
- **D-01 (ii) 重新有守卫**：轮上限提取为具名 `grow::driver::max_rounds(max_loops, max_n_steps)` 并单测（顺带修掉 `+ 1` 不饱和的 debug 溢出）；投降路径本身由 `grow::tests::driver::an_unsatisfiable_hard_core_terminates_and_says_so` 钉住（4 × 5 珠、严格硬核、一遍预算，毫秒级）。
- **退休 `AtomRestraint::periodic_box`（原记于 2026-09-14 条）**：唯一的 `.inp` 生产者传 `[false; 3]`，公开面无生产者，整条 `derive_periodic_box` 路径只有 molpack 自己的测试能走。删除：trait 方法与两处转发、`InsideBoxRestraint::periodic` 字段与两个构造器参数、`cube_from_origin`、`derive_periodic_box`、`PackError::ConflictingPeriodicBoxes` 及 Python 的 `ConflictingPeriodicBoxesError`（wheel 公开面的破坏性变更，stage = experimental）。周期性只在入口声明。`zero_extent_declaration_is_rejected` 改走活的生产者（`with_periodic_box`）。
- **D-05 清**：`initial()` 改收 `&mut RigidView`（调用方本来就握着它），两处 `RigidView::fresh + copy_from_slice` 过桥删除——少两次 6N 分配，也不再存在"放置向量的副本与被重建的坐标可能不一致"这种第二家。

**Why:** 操作者 2026-09-20：按审查结果修复，不留欠债与违反原则的东西。
**How to apply:** 网格覆盖只有 `coverage_radmax` 一个答案；任何"这个 restraint 是不是周期的"问题都问入口，不问形状；生长的格点查询一律经 `CellGrid`。

## 2026-09-20 — molrs 删掉了 shared-dylib 链接形式，molpack 的那一半残留一起清

- **删除**：根 `build.rs` 与 `python/build.rs`（只为动态链接注入 libstd rpath）、`molrs-ffi` 这条 **dev-dependency** 与它唯一的引用点（集成测试 `abi_line`，存在的意义就是把测试二进制拉进「共享 dylib 圈」）、CI 的 `link-dynamic` job（调用 molrs 已删除的 `scripts/verify-shared-dylib.sh`）、两个 manifest 里「七个 native root 必须逐字节一致」的注释与 `sync-dylib-locks.sh` 引用。
- **保留**：wheel 侧 `python/Cargo.toml` 的 `molrs-ffi` 是**真依赖**（capsule 句柄层），ABI 线规则（law P10）与三道 gate 原样不动；`[profile.release]` 两处仍需手工保持一致，只是理由从「dylib 指纹」变成「发布 wheel 与本地 release 同优化」。

**Why:** molrs `76f9bf1d build!: remove the shared-dylib link form and its machinery` —— 那条可选动态链路没有消费者（molpack 一直静态链接），代价是两个仓库里一整套脚本、workflow、feature 与注释。
**How to apply:** 不要再提「链接形式」；molpack 静态链接 molrs，唯一跨进程契约是 capsule 的 major.minor 线。

## 2026-09-20 — one test tier: unit, in-module, owned; no e2e, no regression, no bench

- **`tests/`、`benches/`、`regressions/` 全部删除。** 行为测试只剩 `src/` 里的 `#[cfg(test)]` 模块，与拥有该行为的代码同住（law § 11）；测试体过大就开子模块（`src/grow/tests/{internal,field,prior,entry}.rs`、`src/pipeline/tests.rs`）。原 17 个集成文件 ~236 个 `#[test]` 收敛为 351 个 lib 单测 + 21 个 doc 测试，整层跑完 < 1 s。
- **被删的是什么。** 端到端打包场景（`packer.rs`、`triclinic.rs`、`cli.rs`、`examples_batch.rs`、grow 的 27 个整链/统计用例）、逐位连续性对照（pipeline ≡ preset / `with_restart`）、所有硬编码 golden（`*_regression_*_golden`，Rust 与 Python 两侧）、criterion benches 与 `bench.yml`、`mt_scaling` 测量 harness。判据：断言的是「这次构建碰巧产出的数字」或「整条流水线跑通」，而不是某个模块自己拥有的性质。
- **补的缺口。** 三处原来只被 e2e 覆盖的规则改成属主处的单测：`system::build`（short radius 必须更短，按 target/atom 命名）、`entry::setup`（global restraint 的广播等价律）、`pipeline`（`NoTargets`）。
- **Python 侧同规矩。** 绑定层只测 marshalling / 错误映射 / 回调是否真的回到 Python；`test_integration.py`、`test_pack_peo_examples.py`、`test_pack_peo_topo.py`、`test_fixed_only.py` 删除，其余文件里重复 Rust 物理的用例删除（196 → 143）。

**Why:** 操作者 2026-09-20 裁定：与 molrs 对齐（molrs 已无 tests/ benches/ examples/），且打包质量的度量系统要重新设计，旧的 e2e/回归/bench 既慢又是「守住当前数值」而非守住契约。
**How to apply:** 新行为的测试写进属主模块；想加端到端场景 → 写成 `examples/` 里可跑的程序，不要进测试层；想加 benchmark → 等新度量系统的 spec，不要临时重开 `benches/`。

## 2026-09-14 — regions are molrs's; molpack lifts, never describes

- **One region vocabulary.** A region — sphere, box, triclinic cell, half-space, cylinder, ellipsoid, mesh-bounded solid, union of spheres, and any `&` / `|` / `~` composition — is a `molrs::core::Region` (signed `distance`, negative inside). Every shape describes its inside; outside is `NotRegion`. molpack has **no** `*Region` type, no `StlRegion`, no BVH, no file entry.
- **One lift.** `RegionRestraint(Arc<dyn Region + Send + Sync>)` is the only public geometric restraint (`scale·max(0, distance)²` — distance-quadratic, so `scale` like the `.inp` box kernel); `CellRestraint` is the same lift over a `Parallelepiped` plus the sole `declared_cell` producer. Python accepts any object exposing `_ffi_regionref_capsule()` (`molrs.RegionRef/<line>`, name from `molrs_ffi::abi`, never hard-coded).
- **Packmol-parity kernels are `.inp`-private.** `Inside*/Outside*/Above*/Below*` structs live in `restraint::geometric` (`pub(crate)`), built only by `script::build`; their tests sit beside them in-module. The Python classes of the same names are gone. `AtomRestraint::periodic_box` 已于 2026-09-20 连同 `derive_periodic_box`、`InsideBoxRestraint::periodic` 与 `ConflictingPeriodicBoxes(Error)` 一并删除。

**Why:** the PEO-in-void case needs a region built from atoms in memory (`~SphereUnion`), not a mesh file; and three parallel region models (molrs boolean, molpack SDF + mesh, molvis TS) for one concept violated §1/§9. Operator rulings 2026-09-14: fold into molrs's unreleased 0.14; method name `distance`; one type per shape, sign by `Not`; thorough removal from the Python/Rust public surface.
**How to apply:** a new *shape* is a molrs `impl Region`; a new *penalty* is a molpack `AtomRestraint`; `.inp` grammar untouched (P2). Supersedes specs `stl-region-01/02` and the triclinic DRAFT's `InsideCellRegion`.

## 2026-09-14 — 债务 D-06：公开 region 不受"平面横跨周期轴"拒绝的保护 — 已清（2026-09-20）

`entry::setup::reject_planes_across_periodic_axes` 只读 `AtomRestraint::plane_normal`，而唯一的生产者是 `.inp` 私有的 `Above/BelowPlaneRestraint`；经 `RegionRestraint` 提升的 molrs `HalfSpace`（或任何无界 region）不声明法向，故一个法向沿周期格矢的半空间会被静默接受，其"内侧"随原点位置而变。原三斜集成测试里两条针对公开面的拒绝测试已退役，守卫钉在 `entry::setup::periodic_declaration_tests`。
**Why:** law § 5（隐藏决策、暴露契约）——region 的"无界方向"是 molrs 的几何事实，不应由 molpack 猜；给 `Region` 加"无界轴"查询是一次 molrs 侧的接口决定，超出 region-restraint 链的范围。
**How to apply:** 决定 molrs `Region` 是否暴露该查询后，让 `reject_planes_across_periodic_axes` 同时读它；在此之前，文档已写明公开 region 不受此检查（`docs/rust/restraints-and-pbc.md:81`），所以这是**可见的例外**，不是隐形债。

**2026-09-20 裁定与落地（操作者选 (a)，名字定为 `repeats_along`）**：

- **molrs**：`Region::repeats_along(shift: [F; 3]) -> bool`，默认 `false`（没想过的形状不声称自己重复）。`HalfSpace` 实现为「`shift` 与法向正交」；`And`/`Or` 要求每个成员都重复，`Not` 与被补的区域同真假；有界形状一律走默认。共享契约检查 `check_region_contract` 里加了 `check_repeat_claims`：**凡是声称重复的，每个探针点平移 1x/2x/−3x 后必须还在同一侧**——声称本身带来义务，不声称就不检查。molrs 2068 测试全绿。
- **molpack**：`AtomRestraint::plane_normal()` 删除，换成把规则说出口的 `holds_along(shift) -> bool`（与 molrs 的 `repeats_along` 同形；「hold」是「仍然成立」，不是「关住」——关住只是两个充分条件之一），默认 `true`（看不透的用户 `f`/`fg` 按其他属性一样的规矩，取其自陈）。`RegionRestraint` 用 `repeats_along(shift) || bounded_along(shift)`（后者=AABB 在 shift 触及的轴上有限）回答；`.inp` 的 `Above/BelowPlaneRestraint` 用共享的 `geometric::plane_repeats_along` 回答。`reject_planes_across_periodic_axes` → `reject_restraints_across_periodic_axes`，只读这一个契约；错误 `PlaneAcrossPeriodicAxis { axis, normal }` → `RestraintAcrossPeriodicAxis { axis, restraint }`。
- **为什么是「重复 或 有界」两条**：restraint 在实验室坐标上求值、从不 wrap，所以它必须在每个镜像里给同样的答案。把原子关在一个镜像里（盒/球/cell）算；跟着格子重复（沿面内方向的半空间）也算；两者都不是的才是开口的。副作用是好的：`Cuboid & HalfSpace` 的交集有界，于是「盒子里的一个面」可用；`~Sphere`（void 场景）照旧可用。
- **覆盖**：`entry::setup::region_under_wrap_tests` 五条——提升的半空间横跨周期轴被具名拒绝、同一 region 在受限轴上就是 slab、有界 region 每轴可用、有界成员救回开口成员、void 不受影响。

## 2026-08-28 — solver 设计四原则（用户裁定）

1. **molpack 是纯几何的**：packing 算法（Solver 接缝上的一切）永不依赖力场；构象先验只收用户提供的几何数据（扭转态权重、C∞、持久长度、模板值）。in-loop optimizer（`GenCanPack::with_optimizer`，取代原 relaxer；2026-09-29 起不再受 `ff` 门控）是可选增强，不受此条约束也不得成为 solver 依赖。
2. **packer(GENCAN) 与 grow 平级**：solver 不得调用 pack 内部代码（`pgencan` / `run_phase` / `run_iteration`），但共享同一套架构与生命周期——基础设施段 ①②⑤、`PackContext`、共享 objective、冻结公开 `State`。
3. **用户按 target 选择方法**（`Target::with_method`）：molpack 不替用户判断、不静默退化（小分子不自动降为刚体）；不支持的组合报具名错误。
4. **算法须同时适配 all-atom 与 CG**：排除深度、角度处理、可旋转键感知等决策不得硬编码 AA 假设。

**Why:** 用户在 chain-growth-solver spec 评审后逐条裁定（原话："这个必须牢记" 指第 2 条）；已合并进该 spec 的「设计原则」节。
**How to apply:** 任何触及 `src/stage.rs` / `src/grow/` / `Target::with_method` / packer dispatch 的 spec 或实现改动，先对照四条再动手。待 chain-growth-solver 落地后本条是 CLAUDE.md Hard rules 的晋升候选。

> 2026-09-02：上条的四原则已作为项目法则收进 `.claude/notes/law.md`（ids `pure-geometry-solvers` / `solvers-are-peers` / `user-picks-method` / `aa-and-cg`），CLAUDE.md 各保留一行索引；本条留作历史与理由。
> 2026-09-04：P7 公开名随 `result-as-state` 更新为冻结 `State` 与 `with_restart`（不再写 `PackResult` / `seeded_from`）。`src/solver.rs` 已改名为 `src/stage.rs`。

## 2026-09-02 — packing = 初始构象构建；按不可修复度阶梯分族

- **前提**：键长不必落到平衡值（声明容差，建议 10–15%），后面一定跟力场 minimize；算法只消费键图与几何，化学（元素、键级、反应、氢、立体中心）只能作为用户数据从边界进入。
- **阶梯**：L0 连接性 → L1 拓扑态（链环 / 穿刺 / 贯穿）→ L2 大尺度统计 → L3 密度均匀 → L4 局域重叠 → L5 键长键角。packer 对 L0–L3 负全责（构造保证或硬约束），L4–L5 只承诺在容差内并如实报告残余。
- **五个正交族**（按对状态的操作分）：生成 / 连接 / 精修 / 守卫 / 装饰。"放置"不是族：GenCanPack = 生成(给定构象) + 精修⟨刚体⟩。CBMC 是步选择轴、格是空间轴、硬核是排除体积轴——一个生成器是六轴元组（空间、排除体积、步选择、死路策略、调度、分辨率），先验是数据不是轴。
- **架构**：所有阶段共享 `Stage` trait（requires / guarantees / validate / run），`Pipeline` 拥有生命周期与 handler；族内可替换件各一个 trait（`Space` / `ExcludedVolume` / `Term` / `Pairing` / `Invariant`）；状态 `PackState` 只含拓扑与笛卡尔坐标，刚体自由度是 GENCAN 阶段的局部视图。`Pipeline` 现已 earn（多阶段配方、`Repeat`、`Guarded` 是调用方）——supersede engine-entry-split 的「不引入 Pipeline」Non-goal。
- **熔体主路径**改为理想链生成 + 笛卡尔距离几何精修（Auhl 路线），硬核生长退为受限 / 刷 / 环等需要生成期排除体积的场景。

**Why:** grow 族算法评审（连续版硬核构造在熔体密度下不可行；刚体 GENCAN 做不了 push-off：容差爬坡止步 0.6 Å；生长偏置收缩无修复）+ packing 分类 rev 2，用户逐条裁定。
**How to apply:** specs `stage-pipeline` / `dg-refine` / `grow-axes`（`.claude/specs/`）；任何新 packer 的 spec 须标注它守住的阶梯层，且不得要求逐位精确的成键几何。

## 2026-09-02 — 债务 D-01：连续生长的硬核在熔体密度下不可达；预算耗尽后全局缩核自由落体

- **现象**：当时的集成测试 `grow_cg_kremer_grest_c_inf`（KG 熔体夹具，2026-09-20 随 `tests/` 删除；150 × 100 珠 KG，ρ* = 0.85，tolerance 0.85σ）断言 `fdist == 0 && softened == 0`，实测 fdist = 0.2596（最近对 0.6804σ = 0.8 × 0.85，软化下限）。
- **根因（debugger 只诊断，2026-09-02）**：`src/grow/driver.rs` 的 `hard_scale` 是全局单调不恢复的标量；`regrow_budget = max_loops × n_chains`（:148）只随链数不随链长标度，在第 3365 轮耗尽后 `|| regrow_events >= regrow_budget` 使**任何**死路都缩核 3%，一轮内 0.9127 → 0.8000（8 次 softened 中 5 次在同一轮）；`force_place` 未触发（0 次），提交冲突全部回滚（175 次）。指数回撤 `retract·2^(deadends/4)` 在 40 次死路时达 10240 步，把整条链撤到 stage 0 后在 75% 填充的盒里播种饥饿。
- **反事实**：max_loops = 2000 仍失败（最近对 0.687）；关闭 `soften_after` 只走预算路径 → 8 次缩核集中在一轮；`min_hard_scale = 1.0`（严格硬核）→ 活锁 621 475 轮，稳态 11 630 / 15 000 原子（ρ ≈ 0.66），35 条链永久停在 stage 0；trials = 64 也只到 97.8%。`driver.rs:158` 的轮循环无上限。WLC 角先验 + exclusion_depth 2 不是原因（1-4 自阻塞仅 4.4%/trial，分子间阻塞占 76%）。
- **裁定**：断言本身不可达，属评审 B1（软化按链局部化，归 `grow-axes`）。本地可修的两点先修（`/mol:debug` 最小补丁）：(i) 预算耗尽后的缩核必须走同一条按链 `soften_after` 阶梯，不得每死路一次；(ii) 轮循环加上限，到限即按 abort 契约强制完成并 `converged = false`。预算按步数标度（`max_loops × n_chains × n_steps`）与按链可恢复软化留给 `grow-axes`。测试改为断言算法真正构造性保证的不变量（最近对 ≥ `min_hard_scale × tolerance`、无低于下限的对、缩核为阶梯非自由落体、有限终止），并保留 c_n = 1.76 ± 10% 与角度采样断言；严格硬核保证的那一半标注 TODO(grow-axes ac-004)。

- **已落地（2026-09-02，/mol:debug 最小补丁）**：(i) 缩核只走按链 `soften_after` 阶梯，每轮至多一级（`rung_this_round`）；(ii) 轮循环上限 `max_loops × (max n_steps + 1)`，到限走既有 abort 路径，`force_place` 计入 `softened`；`regrow_budget` 整体删除，`max_loops` 在生长里只剩「每链步数的倍数」一种含义（`src/grow/driver.rs`、`src/stage.rs` 的 `Budget` 文档）。测试当时改为断言真正构造性的不变量（集成测试 `grow_cg_kremer_grest_c_inf` / `grow_dense_strict_core_terminates_unconverged` / `grow_softening_needs_repeated_dead_ends`，已随 `tests/` 于 2026-09-20 删除；现由 `src/grow/tests/driver.rs::an_unsatisfiable_hard_core_terminates_and_says_so` 与 `src/grow/driver.rs::{rung_due_uses_watermark, missed_rung_still_due, force_due_reads_streak_at_floor}` 覆盖）。fast tier 13 目标全绿。**仍欠 `grow-axes`**：按链可恢复软化、预算按步数标度；`src/grow/moves.rs::force_place` 的过期注释已在 2026-09-03 的 docs Mode A 里改为「last-resort placement」措辞，不再欠。

**Why:** 基线 fast tier 有一个红测试挡在 stage-pipeline 之前；法则 §10 禁止跳过或削弱，只允许本地修或路由。
**How to apply:** `grow-axes` 落地时必须吸收本条的 (i)(ii) 与预算标度。**2026-09-20**：KG 熔体夹具随 e2e 删除，但 (i)(ii) 的守卫已在属主处重建——`rung_due` / `retract_depth` / `force_due` / `max_rounds` 各有单测，投降路径由 `grow::tests::driver` 钉住。余下的 `grow-axes`（按链可恢复软化、预算按步数标度）是**待设计的功能轴**，不是欠债；其验收随 spec 一起写。

## 2026-09-02 — 债务 D-02：网格 `radmax` 有两种推导 — 已清（2026-09-20）

- push-off / 生长入口安装网格用 `radmax = max(sys.radius)`（`src/gencan/entry.rs:185`、`src/grow/entry.rs:138`），`initial()` 用 `radmax = 2·max(radius_ini)`（`src/initial.rs:456-461`）；`cell_side = discale·1.01·radmax`（`initial.rs:721-724`），前者的 ±1 模板覆盖约 `1.01·discale·R`，而配对截断是 `2·discale·R`——同一事实两种推导，且前者疑似覆盖不足。
- **裁定**：`stage-pipeline-05` 只保证同一状态下两条拼写用同一种推导（逐位一致），不修正推导本身；修正走 `/mol:debug`（先证明覆盖不足是否真实漏对，再统一到一处）。

**Why:** 架构师在 stage-pipeline 链第三次 design-mode 里发现（law § 9 / § 10）。
**How to apply:** 任何触及 `install_simbox_and_grid` 调用点的改动先看本条；修正落地后删除本条并在 05 的 rustdoc 里去掉引用。
**2026-09-03 更新（05 落地）**：三个阶段前奏的网格 `radmax` 统一改从 `radius_ini` 推导。
**2026-09-20 清账**：覆盖半径本身也统一了——`initial::coverage_radmax`（直径）是唯一答案，半径版确实漏配对，RED 测试见 `initial::grid_coverage_tests`。

## 2026-09-03 — 决定：模板错误的报告顺序随 `Topology` 叶子改变

`stage-pipeline-01` 落地后，`Topology::from_frame` 在读原子数的同时读键表，因此一个**既**少于 3 个原子**又**无键的模板，报错从 `TemplateTooSmall` 变为 `Topology(NoBonds)`；顺序由 `NoAtomsBlock → TemplateTooSmall → NoBonds → RingTemplate` 变为 `NoAtomsBlock → NoBonds → TemplateTooSmall → RingTemplate`。两者都是具名错误，现有测试用例的变体不变；不为保留旧顺序在叶子上加只读原子数的第二个读取器。
**Why:** 一个事实一个家（键图读取只在叶子），叶子不知道 grow 的"至少 3 原子"规则。
**How to apply:** 文档与错误消息以新顺序为准；`src/grow/tests/entry.rs::grow_rejects_two_atom_bondless_as_no_bonds` 钉住该行为（原集成测试层已于 2026-09-20 删除）。

## 2026-09-03 — 债务 D-03：docs/ 的三个 doctest 失败与 rustdoc 断链（预存，非本链引入）— 已清（2026-09-03）

**根因**：`docs/getting_started.md` 的 MkDocs 内容标签（`=== "…"`）体是 4 空格缩进，CommonMark 视为缩进
代码块，rustdoc 当作 Rust 编译；`docs/extending.md` / `docs/concepts.md` 仍写 `Relaxer` /
`RelaxerRunner` / `TorsionMcRelaxer`——这些 trait 已无后继，被 molrs `Optimizer` +
`GenCanPack::with_optimizer` 取代；`src/lib.rs:73-89` / `src/script/` 的链接指向 `ff` / `io` 门控项，
默认特性下不可解析；`pack_context.rs:108` 链接门控模块；`lattice/mod.rs:12` 链接私有 `decorate`。
**处置**：内容标签改为顶层 ```` ```python ```` 围栏；`extending.md` 的 Relaxer 段重写为可编译的
`Optimizer` 示例（`ff` 门控接线段标 `ignore`）；门控项一律改为代码跨度并注明特性；`Relaxer` 用语全站
改为 in-loop optimizer / handler。结果：`cargo doc --no-deps` 0 警告，`cargo test --doc` 16 过 2 忽略。
**余项 — 已清（2026-09-29）**：`docs/rust/handlers-relaxers.md` 改名重写为 `docs/rust/handlers-optimizers.md`，
`docs/rust/index.md`、`zensical.toml`、`docs/index.md` 不再出现 relaxer 用语。

## 2026-09-03 — 债务 D-04：本地 Python 门无法解析环境（molrs 双重 pin）

`uv run --directory python --group typecheck …` 与 `--group dev`（tox）都在解析阶段失败：`molcrafts-molpy` 0.14.0（同级 `../molpy`）把 `molcrafts-molrs` 钉为 `git+https://github.com/MolCrafts/molrs.git@dev#subdirectory=molrs-python`，而 `python/pyproject.toml` 的 `[tool.uv.sources]` 把它钉为 `path = "../../molrs/molrs-python"`，uv 拒绝冲突 URL。chain-growth-solver Task 11 落地记录里已提到同一问题（当时靠 `maturin build` + `uv pip install` 绕过）。
**Why:** 法则 P3——molrs / molpy 的 pin 由人手工管理；harness 不得自动改 pin。它使 `mol_project.build.check` 的 Python 段与 `ci.local` 的 tox 段在本机不可运行，stage-pipeline-05（`python/src/entry.rs` 一行 import）与 -07（Python 镜像）的 Python 验收只能在 CI 或修好 pin 后验证。
**How to apply:** 由 owner 统一 `../molpy` 与 `python/pyproject.toml` 对 molrs 的 pin（同为 path 或同为 git）；在此之前，`/mol:impl` 对 Python 验收项标注"本环境不可运行，待 CI"。
**2026-09-03 更新（07 落地前）**：直接调用 `tox -c python -e py`（prek pre-push 的拼写）**在本地能通过**——tox 环境用 pip
从路径安装 sibling molrs / molpy，pip 不读 molpy 的 `[tool.uv.sources]` git 源，绕开了 uv 的 URL 冲突；
158 个 Python 测试全绿。顺带修了 `python/pyproject.toml` tox `commands_pre[3]` 的优先级 bug
（`Path(..)/'version.py'.read_text()` 先对字符串调用 `.read_text()`，加括号），该 bug 也会让 CI 的
tox 步骤在越过 uv 后失败。余下的债只剩 `uv run --directory python --group dev tox -e py` 这一拼写
（law P4 与 CI 用它）仍因 molpy 的 git pin 无法解析——统一 pin 仍归用户。
**2026-09-14 更新（用户裁定）**：本地开发**没有 pin 这回事**——同一个 venv 里 molrs / molpy / molpack 三个包都是 sibling 源码树的 editable 安装，用 `--no-deps` 装，不做任何跨包解析（`maturin develop --release` 装 molrs；`uv pip install --no-deps -e ../molpy`；`uv pip install --no-deps --no-build-isolation -e python/` 装 molpack，之后 `maturin develop --release --skip-install` 原地重建 `.so`）。`[tool.uv.sources]` 只服务 CI / 发布。`python/.venv` 曾装着 9 月 10 日的非 editable 拷贝，遮住了新构建；已重装为 editable。**不得**再把"pin 冲突"当本地阻塞项报告，也不得用 PYTHONPATH 绕。


## 2026-09-03 — 发现：`PlacementsMut` 的 COM 访问器无越界保护

`set_com(nmol, ..)` 的偏移 `3*nmol` 正好是分子 0 的 Euler 槽，不 panic 而静默串块；`set_euler(nmol, ..)` 才越界 panic。`stage-pipeline-02` 的 `RigidView` 访问器加显式 `i < nmol` 断言并有测试钉住。
**Why:** law § 8（非法状态不可表示）；tester 在写 02 的 RED 时用参考桩实测发现。
**How to apply:** 任何按 (COM 块 | Euler 块) 布局索引的代码都要显式检查分子下标，不能依赖偏移算术。

## 2026-09-03 — 债务 D-05：`initial()` 仍讲扁平切片，借临时 `RigidView` 过桥两次 — 已清（2026-09-20）

**现象**：02 把 xcart 重建收进 `RigidView::write_xcart` 后，`src/initial.rs` 的两处内部重建
（首猜后、Phase-1 后）各自 `RigidView::fresh(ntotmol)` → `copy_from_slice(x)` → `write_xcart`，
因为 `initial(x: &mut [F], ..)` 的 460 行主体全是 `ilubar / ilugan / icart` 偏移算术，唯一真实
调用方 `gencan/solver.rs::GencanStage::run`（当时名 `GencanSolver::solve`）手里其实握着 `&mut RigidView`。
**为何不在 02 改**：改签名只是把摩擦搬成三处 `let x = view.as_mut_slice();` 重绑，主体不变；
两次 6N 堆分配不在热路径。02 的 Files / Tasks 不含 `initial()` 的重构（law § 4）。
**归属**：05——`GencanStage::run` 前奏接管 `prepare()` 后，`initial()` 或改收 `&mut RigidView`
整段持有一个视图，或被前奏吸收；届时删除这两处过桥。
**2026-09-20 清账**：`initial(view: &mut RigidView, ..)` 落地，两处过桥删除，体内按需 `let x = view.as_mut_slice()` 重绑。

<!-- mol:note:topic:result-as-state -->
## 2026-09-04 — 公开冻结结果叫 State；接续只留 with_restart

`PackResult` 就地改名为 `State`（Rust crate + Python pyclass）。直播 `PackState` 不改名、不与 `State` 合并。`GenCanPack::seeded_from` 删除；跨入口接续只留 `with_restart`。高分子熔体评测主张 ρ = 1.2 g/cm³，不是 `PackEngine::with_density` 的算法上限。
**Why:** 用户裁定一套公开结果；P7 示例必须跟现行公开名，否则法则教已删 API。
**How to apply:** 新代码与文档只写 `State` / `with_restart`。禁止 `PackResult` 别名、禁止公开 `GrowState`。
**Supersedes:** law P7 正文里的 `PackResult` / `seeded_from` 拼写。
