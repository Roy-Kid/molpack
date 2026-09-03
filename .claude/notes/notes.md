# .claude/notes/notes.md — molpack evolving decisions

Use `/mol:note <decision>` to add or supersede an entry. Rules that become absolute move to `law.md` (one line each in CLAUDE.md); this file keeps the reasoning and the history.
Format per entry:

```
## YYYY-MM-DD — short title
<rule>
**Why:** <reason — incident, deadline, or stakeholder ask>
**How to apply:** <when / where this kicks in>
```

## 2026-08-28 — solver 设计四原则（用户裁定）

1. **molpack 是纯几何的**：packing 算法（Solver 接缝上的一切）永不依赖力场；构象先验只收用户提供的几何数据（扭转态权重、C∞、持久长度、模板值）。`ff` 门控的 in-loop relaxer 是可选增强，不受此条约束也不得成为 solver 依赖。
2. **packer(GENCAN) 与 grow 平级**：solver 不得调用 pack 内部代码（`pgencan` / `run_phase` / `run_iteration`），但共享同一套架构与生命周期——基础设施段 ①②⑤、`PackContext`、共享 objective、`PackResult`。
3. **用户按 target 选择方法**（`Target::with_method`）：molpack 不替用户判断、不静默退化（小分子不自动降为刚体）；不支持的组合报具名错误。
4. **算法须同时适配 all-atom 与 CG**：排除深度、角度处理、可旋转键感知等决策不得硬编码 AA 假设。

**Why:** 用户在 chain-growth-solver spec 评审后逐条裁定（原话："这个必须牢记" 指第 2 条）；已合并进该 spec 的「设计原则」节。
**How to apply:** 任何触及 `src/solver.rs` / `src/grow/` / `Target::with_method` / packer dispatch 的 spec 或实现改动，先对照四条再动手。待 chain-growth-solver 落地后本条是 CLAUDE.md Hard rules 的晋升候选。

> 2026-09-02：上条的四原则已作为项目法则收进 `.claude/notes/law.md`（ids `pure-geometry-solvers` / `solvers-are-peers` / `user-picks-method` / `aa-and-cg`），CLAUDE.md 各保留一行索引；本条留作历史与理由。

## 2026-09-02 — packing = 初始构象构建；按不可修复度阶梯分族

- **前提**：键长不必落到平衡值（声明容差，建议 10–15%），后面一定跟力场 minimize；算法只消费键图与几何，化学（元素、键级、反应、氢、立体中心）只能作为用户数据从边界进入。
- **阶梯**：L0 连接性 → L1 拓扑态（链环 / 穿刺 / 贯穿）→ L2 大尺度统计 → L3 密度均匀 → L4 局域重叠 → L5 键长键角。packer 对 L0–L3 负全责（构造保证或硬约束），L4–L5 只承诺在容差内并如实报告残余。
- **五个正交族**（按对状态的操作分）：生成 / 连接 / 精修 / 守卫 / 装饰。"放置"不是族：GenCanPack = 生成(给定构象) + 精修⟨刚体⟩。CBMC 是步选择轴、格是空间轴、硬核是排除体积轴——一个生成器是六轴元组（空间、排除体积、步选择、死路策略、调度、分辨率），先验是数据不是轴。
- **架构**：所有阶段共享 `Stage` trait（requires / guarantees / validate / run），`Pipeline` 拥有生命周期与 handler；族内可替换件各一个 trait（`Space` / `ExcludedVolume` / `Term` / `Pairing` / `Invariant`）；状态 `PackState` 只含拓扑与笛卡尔坐标，刚体自由度是 GENCAN 阶段的局部视图。`Pipeline` 现已 earn（多阶段配方、`Repeat`、`Guarded` 是调用方）——supersede engine-entry-split 的「不引入 Pipeline」Non-goal。
- **熔体主路径**改为理想链生成 + 笛卡尔距离几何精修（Auhl 路线），硬核生长退为受限 / 刷 / 环等需要生成期排除体积的场景。

**Why:** grow 族算法评审（连续版硬核构造在熔体密度下不可行；刚体 GENCAN 做不了 push-off：容差爬坡止步 0.6 Å；生长偏置收缩无修复）+ packing 分类 rev 2，用户逐条裁定。
**How to apply:** specs `stage-pipeline` / `dg-refine` / `grow-axes`（`.claude/specs/`）；任何新 packer 的 spec 须标注它守住的阶梯层，且不得要求逐位精确的成键几何。

## 2026-09-02 — 债务 D-01：连续生长的硬核在熔体密度下不可达；预算耗尽后全局缩核自由落体

- **现象**：`tests/grow.rs::grow_cg_kremer_grest_c_inf`（150 × 100 珠 KG，ρ* = 0.85，tolerance 0.85σ）断言 `fdist == 0 && softened == 0`，实测 fdist = 0.2596（最近对 0.6804σ = 0.8 × 0.85，软化下限）。
- **根因（debugger 只诊断，2026-09-02）**：`src/grow/driver.rs` 的 `hard_scale` 是全局单调不恢复的标量；`regrow_budget = max_loops × n_chains`（:148）只随链数不随链长标度，在第 3365 轮耗尽后 `|| regrow_events >= regrow_budget` 使**任何**死路都缩核 3%，一轮内 0.9127 → 0.8000（8 次 softened 中 5 次在同一轮）；`force_place` 未触发（0 次），提交冲突全部回滚（175 次）。指数回撤 `retract·2^(deadends/4)` 在 40 次死路时达 10240 步，把整条链撤到 stage 0 后在 75% 填充的盒里播种饥饿。
- **反事实**：max_loops = 2000 仍失败（最近对 0.687）；关闭 `soften_after` 只走预算路径 → 8 次缩核集中在一轮；`min_hard_scale = 1.0`（严格硬核）→ 活锁 621 475 轮，稳态 11 630 / 15 000 原子（ρ ≈ 0.66），35 条链永久停在 stage 0；trials = 64 也只到 97.8%。`driver.rs:158` 的轮循环无上限。WLC 角先验 + exclusion_depth 2 不是原因（1-4 自阻塞仅 4.4%/trial，分子间阻塞占 76%）。
- **裁定**：断言本身不可达，属评审 B1（软化按链局部化，归 `grow-axes`）。本地可修的两点先修（`/mol:debug` 最小补丁）：(i) 预算耗尽后的缩核必须走同一条按链 `soften_after` 阶梯，不得每死路一次；(ii) 轮循环加上限，到限即按 abort 契约强制完成并 `converged = false`。预算按步数标度（`max_loops × n_chains × n_steps`）与按链可恢复软化留给 `grow-axes`。测试改为断言算法真正构造性保证的不变量（最近对 ≥ `min_hard_scale × tolerance`、无低于下限的对、缩核为阶梯非自由落体、有限终止），并保留 c_n = 1.76 ± 10% 与角度采样断言；严格硬核保证的那一半标注 TODO(grow-axes ac-004)。

- **已落地（2026-09-02，/mol:debug 最小补丁）**：(i) 缩核只走按链 `soften_after` 阶梯，每轮至多一级（`rung_this_round`）；(ii) 轮循环上限 `max_loops × (max n_steps + 1)`，到限走既有 abort 路径，`force_place` 计入 `softened`；`regrow_budget` 整体删除，`max_loops` 在生长里只剩「每链步数的倍数」一种含义（`src/grow/driver.rs`、`src/solver.rs` 文档）。测试改为断言真正构造性的不变量（`tests/grow.rs`：`grow_cg_kremer_grest_c_inf` 保留 c_n 与角度断言 + 有效尺度上的最近对/阶梯守卫；新增 `grow_dense_strict_core_terminates_unconverged`、`grow_softening_needs_repeated_dead_ends`）。fast tier 13 目标全绿。**仍欠 `grow-axes`**：按链可恢复软化、预算按步数标度；`src/grow/moves.rs::force_place` 的过期注释已在 2026-09-03 的 docs Mode A 里改为「last-resort placement」措辞，不再欠。

**Why:** 基线 fast tier 有一个红测试挡在 stage-pipeline 之前；法则 §10 禁止跳过或削弱，只允许本地修或路由。
**How to apply:** `grow-axes` 落地时必须吸收本条的 (i)(ii) 与预算标度；任何触及 `driver.rs` 死路梯的改动先跑 KG 夹具的 radscale 轨迹测试。

## 2026-09-02 — 债务 D-02：网格 `radmax` 有两种推导

- push-off / 生长入口安装网格用 `radmax = max(sys.radius)`（`src/gencan/entry.rs:185`、`src/grow/entry.rs:138`），`initial()` 用 `radmax = 2·max(radius_ini)`（`src/initial.rs:456-461`）；`cell_side = discale·1.01·radmax`（`initial.rs:721-724`），前者的 ±1 模板覆盖约 `1.01·discale·R`，而配对截断是 `2·discale·R`——同一事实两种推导，且前者疑似覆盖不足。
- **裁定**：`stage-pipeline-05` 只保证同一状态下两条拼写用同一种推导（逐位一致），不修正推导本身；修正走 `/mol:debug`（先证明覆盖不足是否真实漏对，再统一到一处）。

**Why:** 架构师在 stage-pipeline 链第三次 design-mode 里发现（law § 9 / § 10）。
**How to apply:** 任何触及 `install_simbox_and_grid` 调用点的改动先看本条；修正落地后删除本条并在 05 的 rustdoc 里去掉引用。

## 2026-09-03 — 决定：模板错误的报告顺序随 `Topology` 叶子改变

`stage-pipeline-01` 落地后，`Topology::from_frame` 在读原子数的同时读键表，因此一个**既**少于 3 个原子**又**无键的模板，报错从 `TemplateTooSmall` 变为 `Topology(NoBonds)`；顺序由 `NoAtomsBlock → TemplateTooSmall → NoBonds → RingTemplate` 变为 `NoAtomsBlock → NoBonds → TemplateTooSmall → RingTemplate`。两者都是具名错误，现有测试用例的变体不变；不为保留旧顺序在叶子上加只读原子数的第二个读取器。
**Why:** 一个事实一个家（键图读取只在叶子），叶子不知道 grow 的"至少 3 原子"规则。
**How to apply:** 文档与错误消息以新顺序为准；`tests/topology.rs::topology_error_precedence_no_bonds_before_too_small` 钉住叶子侧行为。

## 2026-09-03 — 债务 D-03：docs/ 的三个 doctest 失败与 rustdoc 断链（预存，非本链引入）

`cargo test -p molcrafts-molpack --doc` 失败 3/19：`docs/getting_started.md:20,28` 是未标注语言的 Python 代码块被 rustdoc 当作 Rust 编译；`docs/extending.md:342` 引用已不存在的 `molpack::Relaxer` / `RelaxerRunner` 与 `rand::RngCore`。`cargo doc --no-deps` 另有约 18 条 intra-doc 断链（`Optimizer`、`TorsionMcOptimizer`、`Script::build`、`crate::Relaxer` …）与一条 `lattice` 链接到私有 `decorate` 的警告。全部属于工作树里在途的 `Relaxer → Handler` 改名与 docs 重写，不由 stage-pipeline 引入。
**Why:** law § 10——看见即记录；stage-pipeline 各子 spec 的 docs 门是 `cargo doc --no-deps` 零警告，这条债务不清则 04/05 的 ac 无法按字面通过。
**How to apply:** 路由 `/mol:docs`（Mode A）修 fence 语言标注与死引用；在 stage-pipeline-04 之前清掉，否则其 docs 类验收只能以"新增符号零警告"为准并把预存警告列为例外。

## 2026-09-03 — 债务 D-04：本地 Python 门无法解析环境（molrs 双重 pin）

`uv run --directory python --group typecheck …` 与 `--group dev`（tox）都在解析阶段失败：`molcrafts-molpy` 0.14.0（同级 `../molpy`）把 `molcrafts-molrs` 钉为 `git+https://github.com/MolCrafts/molrs.git@dev#subdirectory=molrs-python`，而 `python/pyproject.toml` 的 `[tool.uv.sources]` 把它钉为 `path = "../../molrs/molrs-python"`，uv 拒绝冲突 URL。chain-growth-solver Task 11 落地记录里已提到同一问题（当时靠 `maturin build` + `uv pip install` 绕过）。
**Why:** 法则 P3——molrs / molpy 的 pin 由人手工管理；harness 不得自动改 pin。它使 `mol_project.build.check` 的 Python 段与 `ci.local` 的 tox 段在本机不可运行，stage-pipeline-05（`python/src/entry.rs` 一行 import）与 -07（Python 镜像）的 Python 验收只能在 CI 或修好 pin 后验证。
**How to apply:** 由 owner 统一 `../molpy` 与 `python/pyproject.toml` 对 molrs 的 pin（同为 path 或同为 git）；在此之前，`/mol:impl` 对 Python 验收项标注"本环境不可运行，待 CI"。
