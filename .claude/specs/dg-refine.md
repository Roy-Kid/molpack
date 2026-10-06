---
title: dg-refine — 笛卡尔距离几何精修阶段（软核 + 键距弹簧 + 半径阶梯）
status: draft
created: 2026-09-02
chain: packing-taxonomy（02 of 3；依赖 stage-pipeline 的 PackState / Stage）
---

# dg-refine

状态：DRAFT（2026-09-02）。依赖 `stage-pipeline`（`PackState`、`Stage`、`Pipeline`）。
与 `grow-axes` 并行，无代码耦合。

**设计来源**：grow 族算法评审 A1 / A2（刚体 GENCAN 在熔体密度下做不了 push-off：
`lattice_ramp` 基准里容差爬坡止步于 0.6 Å）与 packing 分类 rev 2「精修族 R1」。
用户裁定：键长不必精确、后接力场 minimize、全部几何无力场。

## Goal

新增精修族的第一个阶段 `DgRefine`（`Refine<Cartesian, Spg>`）：变量是**全部自由原子的
笛卡尔坐标**；目标函数是可加的几何项之和——既有的分子间重叠项（`objective.rs::pair_term`，
在半径阶梯的当前 scale 上）、成键距离弹簧（1-2 与 1-3 对拉向模板距离，可选 1-4）、可选
签名体积项（用户给定的四元组 + 符号，保手性）、既有逐原子约束项；优化器复用
`gencan::gencan`（SPG + CG 的有界最小化，`src/gencan/mod.rs:130`），梯度复用 objective 已
算出的逐原子 `gxcar`（跳过 `project_cartesian_gradient`，`src/objective.rs:1268`）；调度是
半径阶梯：scale 从 `s0` 线性爬到 1.0，每级收敛到 precision 或用完每级预算再升。

出口保证：L4（重叠到声明容差，尽力而为并如实报告）、L5（键长 / 1-3 距离在
`bond_tolerance` 内）；L0–L2 逐位不变（不改键图；大尺度统计变化在验收里量化并断言 ≤ 阈值）。

可观测行为：熔体密度下"生成 → DgRefine"在 molpack 内部把最近对推到容差，不必出境到 MD；
`State::fdist` 由共享 objective 在 scale = 1.0 上报告；残余如实给出。

```rust
let refined = Pipeline::new()
    .stage(CbmcGrow::new(prior).with_initial_hard_scale(0.5))
    .stage(DgRefine::new().with_ladder(0.5, 10).with_bond_stiffness(1.0))
    .with_density(1.06).with_seed(42)
    .run(&[peo], 60)?;
assert!(refined.fdist < 0.01 || refined.fdist < grown.fdist * 0.1);
```

## Non-goals

- 不做扭转空间变量、不做 NeRF 雅可比：笛卡尔变量对环与网络一视同仁，是本 spec 的要点。
- 不做 MD、不做热浴：零温阶梯最小化 + 可选的哈希流随机扰动（几何退火）。
- 不感知手性、不感知氢：签名体积四元组与"装饰原子最后挂载"都是数据 / 后续 spec
  （`topology-model`）。
- 不做格上 MC 修复（R3）、不做拓扑守卫的键段交叉检查（守卫族，`ring-closure` spec 引入
  `Invariant` 实现；本 spec 只在每级之间留 `Guarded` 挂点）。
- 不替换 `GenCanPack` 的刚体 push-off：R2 仍是刚性构象体系的正确工具。
- 不做 rayon 新并行：pair kernel 已有的并行路径原样复用。

## Public surface

**Rust 新增**

- `src/refine/config.rs`（**叶子**，同 `grow/config.rs` 纪律：只 import molrs）：
  `pub struct RefineConfig { ladder: (s0, rungs), loops_per_rung, k_bond, k_angle, k_dihedral,
  k_intra, chiral: Vec<([usize; 4], Sign)>, anneal }` 与 builders。
- `src/refine/mod.rs`：`pub struct DgRefine`（`PackEngine` + `StageFactory` 预设）与
  builders：`with_ladder(s0: F, rungs: usize)`（默认 `(0.5, 10)`）、
  `with_loops_per_rung(n)`（默认 = `max_loops / rungs`，至少 1）、
  `with_bond_stiffness(k_b: F)`（默认 1.0）、`with_angle_stiffness(k_a: F)`（1-3，默认 0.5）、
  `with_dihedral_stiffness(k_d: F)`（1-4，默认 0.0 = 关）、
  `with_intra_overlap(k_i: F)`（默认 `1.0`；排除表读 `Target.special_bonds`）、
  `with_chiral(quads: Vec<([usize; 4], Sign)>)`（按模板原子索引，逐拷贝广播）、
  `with_anneal(amplitude: F)`（每级开始的哈希流随机扰动，默认 0.0 = 关）。
- `src/refine/terms.rs`：`pub trait Term: Send + Sync { fn name(&self) -> &'static str;
  fn eval(&self, state: &PackState, mode: EvalMode, grad: Option<&mut [[F; 3]]>) -> F; }`
  与内置项（各 ≤ 400 行，超出则一项一文件）：
  - `OverlapTerm`：既有分子间 pair term（`pair_term` 跳过同分子对，`objective.rs:246`），
    在当前 scale 上；
  - `IntraOverlapTerm`：**同分子**、排除深度之外的对，用同一半径与 scale，自带排除表感知的
    cell 遍历（不改 `pair_term` 热循环——`examples_batch` 逐位风险）；没有它，链只有弹簧
    没有分子内斥力，会自穿；
  - `BondTerm`、`Angle13Term`、`Dihedral14Term`：距离弹簧，成键对来自 `molrs::Topology`
    （直接消费，不从 `grow::internal` 取、不 import `topology_for_growth`——避免 refine → grow 的跨族边）；
  - `ChiralTerm`（签名体积，数据）、`RestraintTerm`（既有 `AtomRestraint::f/fg`）。
  所有项无量纲化到 overlap 项的自然尺度（`tol⁴`）：`E_bond = k_b · tol⁴ · Σ ((d − d₀)/d₀)²`，
  1-3 / 1-4 同形；`E_chiral = k_c · tol⁴ · Σ max(0, −sgn·V/V₀)²`。刚度是几何参数，不是力常数。
- `src/refine/cartesian.rs`：`pub(crate) struct CartesianObjective<'a>` 实现
  `objective::Objective`（`evaluate(x, mode, grad)`；`bounds` 用 trait 默认的无界实现——
  `PackContext::bounds` 假定 COM/Euler 布局，`objective.rs:1488-1511`，不可复用），
  `x` = 自由原子坐标平铺（3·N_free）；写入 `ctx.xcart` 自由槽位、重建 cell list、累加各 `Term`。
- `src/refine/minimize.rs`：一个小的 SPG 外循环（BB 步长 + 非单调线搜索），线搜索直接调用
  `gencan::spg::spgls`（`src/gencan/spg.rs:44`，签名已是 `&mut dyn Objective`）；可选
  `gencan::cg::cg_solve` 做牛顿式加速。**不**动 `gencan/mod.rs`（987 行，已超预算），
  **不**调用 `pgencan / gencan()`（其外层含 Packmol 专属逻辑，且 `GencanParams` 绑定刚体语义）。
- `src/refine/ladder.rs`：`pub struct Ladder { s0, rungs }` 与 `Ladder::scales()` 迭代器；
  每级通过 `PackContext::set_radius`（`gencan/phases.rs:186-191` 的同一机制）设
  `radius = s · radius_ini` 并同步 `atom_props`。
- **前置：`src/objective.rs` 行为保持拆分**（1576 行，已超预算；本 spec 要用的
  `accumulate_pair_fg / accumulate_collective_fg / accumulate_constraint_values_and_gradients_from_xcart`
  都是私有，cell 插入融合在 `expand_molecules` Phase B）：拆为 `src/objective/{mod.rs（pair kernel）,
  expand.rs（expand_molecules + project_cartesian_gradient，刚体专属）, parallel.rs（rayon 变体
  ~:1003-1160）, cartesian.rs（新：`pub(crate) fn bin_xcart_into_cells` +
  `pub(crate) fn accumulate_cartesian_fg`，结果留在 `gxcar`，**不**投影）}`。与 DRAFT 的
  `pair-loop-context-split` 触及同一读集——两者合并为一次 objective 重组，先于本 spec 落地，
  `examples_batch` 逐位守门。

**Python**：`DgRefine` pyclass（构造器 + 上述 builders；`with_chiral` 接受
`list[tuple[tuple[int,int,int,int], int]]`）；`.pyi` / `_protocols.py` 同步。**CLI**：无。

## Module placement

架构师 review-mode 意见（2026-09-02）已合并：新目录 `src/refine/`，是 `gencan/`、`grow/`
的**同级**，不在任一之下。

| 文件 | 内容 | 预算 |
|---|---|---|
| `src/refine/config.rs`（叶子） | `RefineConfig`：阶梯、刚度、排除深度、手性四元组、退火 | ≤ 150 |
| `src/refine/mod.rs` | 门面 + `DgRefine` 预设（超 200 行则拆 `entry.rs`）+ `Stage` impl | ≤ 300 |
| `src/refine/terms.rs` | `Term` trait + 内置项（超 400 则一项一文件 `terms/*.rs`） | ≤ 400 |
| `src/refine/cartesian.rs` | 自由原子 ↔ 平铺向量映射、`Objective` impl、cell list 重建 | ≤ 300 |
| `src/refine/minimize.rs` | SPG 外循环，调用 `gencan::spg::spgls`（+ 可选 `cg_solve`） | ≤ 200 |
| `src/refine/ladder.rs` | 半径调度，走 `PackContext::set_radius` | ≤ 150 |
| `src/objective/cartesian.rs`（前置拆分产物） | `bin_xcart_into_cells`、`accumulate_cartesian_fg` | ≤ 200 |

- **依赖方向**：`refine/*` → `stage` / `context::pack_state` / `molrs::Topology` / `objective::cartesian` /
  `gencan::{spg, cg}` / `restraint`。允许 `refine → gencan::{spg, cg}`：原则 2 与
  `chain-growth-solver.acceptance.md:30` 点名禁止的是 `pgencan / run_phase / run_iteration`
  这些刚体驱动，不是线搜索原语；本 spec 的验收 grep 把同一条禁令扩到 `src/refine/`，并
  追加 `gencan::solver / phases / entry / mod::gencan(` 四项。
- **不**从 `grow::internal` 取键图（跨族边）；refine 直接消费 `molrs::Topology`，永不 import `topology_for_growth`。
- `gencan/mod.rs`（987 行）本 spec 一行不加；若将来要复用其外循环，先把
  `tn_linesearch`（:611-987）拆到 `gencan/tnls.rs`，那是独立的 hygiene 改动。
- `Term` trait 住 `refine/terms.rs`，不进 `restraint/`：约束是逐点谓词，项是全局可加能量。
  `Term` 与 `CartesianObjective` Rust-only；Python 只镜像 `DgRefine` 预设。
- 命名：`DgRefine` 是 `Stage`，**不是** `optimizer/` 里的 in-loop optimizer（`GenCanPack::with_optimizer`）——模块文档写明；
  无 "packmol"；泛型标记若将来引入用 `CartesianVars / SpgMinimizer` 而非 `Cartesian / Spg`
  （后者与模块 `gencan::spg` 混视）。
- 无 `ff`：`grep -rn 'molrs::ff\|molrs::optimize' src/refine/` 无命中（crate 内已无 `cfg(feature = "ff")`，`ff` 只透传 `molrs/ff`）。

## Numerical contract

- **梯度正确**：每个 `Term` 的解析梯度与中心差分（h = 1e-5 Å）相对误差 < 1e-6（含 PBC
  最小镜像跨界的成键对、含 fixed 原子的 pair）；`CartesianObjective` 总梯度同样通过。
- **不改拓扑与统计**：键图逐位不变；对生成阶段产物做全阶梯精修后，每链 R_g 相对变化
  ≤ 3%，内距曲线 ⟨R²(s)⟩/s 在 s ≥ 50 相对变化 ≤ 5%（L2 守恒是本阶段的核心承诺）。
- **成键几何**：精修后 1-2 距离相对模板偏差 ≤ `bond_tolerance`（默认 0.10），1-3 ≤ 1.5×。
- **分子内接触**：排除深度之外的同分子对在末级最小距离 ≥ tolerance − 1e-9 的比例 ≥ 99%
  （独立复算，不走 `fdist`——`fdist` 按 Packmol 语义只计分子间对，本 spec 不改这把尺）。
- **裁决**：末级在 scale = 1.0 上由共享 objective 评估 `fdist / frest`；本阶段不自报。
- **确定性**：同 seed 逐位一致；退火扰动走 `stream(seed, atom, rung, 0, SALT_ANNEAL)`。
- **单调性不作保证**：级间 fdist 可能上升（半径变大）；记录每级 fdist 到 `StepInfo`
  （`loop_idx` = 级号，`radscale` = 当前 scale）。

## Test plan

- `DgRefine` 属主模块内的 `#[cfg(test)] mod tests`（新，default feature；crate 无 `tests/` 目录，共享夹具走 `src/testutil.rs`）：
  1. 各 `Term` 有限差分梯度（含跨周期边界的键、含 fixed 邻居）；
  2. 扰动链恢复：对模板链加 0.3 Å 随机扰动，仅 Bond + Angle13 项精修后成键几何回到
     容差内；
  3. 重叠消除：8 × 12 珠链在 26 Å 盒里人工重叠放置后
     `DgRefine` 到 `fdist < precision`；
  4. **刚体做不到、精修做得到**：重建原 `lattice_grow_then_seeded_push_off_dense`（已随 `tests/` 于 2026-09-20 删除）的
     20 × 24 珠 / 22 Å 算例，断言 `DgRefine` 末级 `fdist ≤ 0.1 × grown.fdist` 且严格小于
     `GenCanPack::with_restart` 同预算的结果；
  5. L2 守恒：精修前后每链 R_g 相对变化 ≤ 3%；
  6. 约束项：`InsideSphereRestraint` 下精修不把原子推出球（`frest == 0`）；
  7. 手性：给定四元组符号，精修后签名体积符号不变；
  8. 确定性与 `should_stop` 中途终止；
  9. 分子内：两条人工折叠到自穿的链（1-6 对 0.5 Å）经 `IntraOverlapTerm` 精修后，各自
     `Target.special_bonds` 表外的同分子最小距离 ≥ tolerance（两张 Target 各一张表：
     深度 1 与深度 3；原则 4）。
- `examples/pack_peo`（io）增加 `refine` 模式：`CbmcGrow(s₀=0.5) → DgRefine` 在
  ρ = 1.06、dp = 100 × 25 上报告每级 fdist、末级最近对、R_g 变化、耗时——scientific
  验收的数据来源；不进 fast tier。
- 五个 `cargo run --release --features io --example pack_<name>` 程序收敛（原 `examples_batch` harness 已于 2026-09-20 删除）（objective.rs 抽取是行为保持搬移）。
- `python/tests/test_refine.py`：builder 链、两阶段管线冒烟、`with_chiral` 类型检查。

## Doc plan

- `CLAUDE.md` 架构表加 `src/refine/` 行；回归命令不变。
- `docs/python/guide/growth.md`「Reading softened」与「Melt density」两节改写：push-off
  的推荐做法是 `DgRefine`，`GenCanPack::with_restart` 保留为刚性构象体系的工具；
  新节「Refining with distance geometry」（阶梯、刚度、手性数据、残余的含义）。
- `docs/python/api-reference.md`、`docs/rust/handlers-optimizers.md`（`StepInfo` 在精修阶段
  的字段语义）。
- `src/refine/mod.rs` 模块文档：与 Crippen–Havel / ETKDG 的关系、与 Auhl push-off 的关系、
  "无力场"边界（刚度是相对 overlap 项的无量纲几何参数）。

## Risks / open questions

1. **零温阶梯的残余**：最后几级可能停在局部极小；缓解是更多级数与退火扰动。PEO
   ρ = 1.06 的残余分布是本 spec 的科学验收项，阈值在 `examples/pack_peo refine` 跑出
   数据后定（初值：末级最近对 ≥ 0.8 × tol 的原子对占比 ≥ 99%）。
2. **弹簧刚度默认值**：`k_b = 1.0`、`k_a = 0.5` 是占位，按 PEO 与 KG 两个算例标定
   （L5 在容差内且 L2 守恒两条同时满足的最小刚度）。
3. **cell list 重建成本**：每次 evaluate 重建 cell list（今天 GENCAN 靠 `x` 键控缓存）；
   精修的每次评估都动全部原子，缓存无用，重建是 O(N)，可接受；若基准显示占比过高，
   改为按位移阈值的增量重分箱。
4. **fixed 原子**：变量只含自由原子；fixed 原子进 pair 项但不进变量；`Topology` 的键图
   含 fixed 拷贝的键但弹簧只对至少一端自由的键生效。
5. **与 `pair-loop-context-split` 的顺序耦合**：objective.rs 的拆分是本 spec 的前置，与该
   spec 触及同一读集，合成一次重组先落地（架构师裁定），再做本 spec。
6. **对接缝的依赖**：本 spec 写在"接缝 trait"上，不硬编码 `Stage`；若 stage-pipeline 延后，
   `DgRefine` 先以 `Solver` + `PackEngine` 入口形态落地（同 `CbmcGrow` 今天的形态），
   接缝升级时只改类型名。
7. **纯最小化的自穿**：`IntraOverlapTerm` 只能推开、不能解开已穿过的链段（L1 不变量），
   这与设计一致——生成阶段负责不制造 L1 缺陷，精修只处理 L4–L5。
