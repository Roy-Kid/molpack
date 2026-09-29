---
title: grow-axes — 一个生长驱动，六条正交轴；CbmcGrow / LatticeGrow / WalkGrow 成为预设
status: draft
created: 2026-09-02
chain: packing-taxonomy（03 of 3；依赖 stage-pipeline 的 PackState / Stage）
---

# grow-axes

状态：DRAFT（2026-09-02）。依赖 `stage-pipeline`。与 `dg-refine` 并行，无代码耦合
（`WalkGrow` + `DgRefine` 的端到端熔体验收在 dg-refine 的 Test plan 里）。

**设计来源**：packing 分类 rev 2「生成族的正交轴」（用户裁定：LatticeGrow 也含 CBMC
思想，方法必须正交）；grow 族算法评审 A3（relax 死功能）、B1（全局软化）、B2（relax
贪心判据）、B3（格相全局 RNG、顺序生长）、C1（`serial` 的 any() 折叠）。

## Goal

把两套生长循环——`src/grow/driver.rs` + `moves.rs`（连续）与 `src/grow/lattice/saw.rs`
（金刚石格）——合成**一个**泛型驱动 `Grow<S: Space, X: ExcludedVolume<S>>`，六条轴各自
独立可换：

| 轴 | v1 取值 | 形态 |
|---|---|---|
| 空间 `Space` | `Continuum`、`Diamond` | trait：给定参考原子与先验抽样产生候选（连续：NeRF 位置；格：3 个非回头方向 + RIS 权重）；`to_continuum` |
| 排除体积 `ExcludedVolume<S>` | `OverlapField`（硬核 + 软壳）、`Occupancy`（格点 + 近邻守卫）、`NoExclusion` | trait：`probe / insert / remove / nearest` |
| 步选择 `Selector` | `Rosenbluth { trials: Sample(k) \| Enumerate, beta }` | enum；候选权重 = 先验权重 × exp(−β·penalty) × [未被阻塞]；`Sample(1), β = 0` 即理想行走，`Enumerate, β → ∞` 即今天的格上走法 |
| 死路策略 `Escape` | 有序梯：`Retract { depth }`、`Reseed { void_bias }`、`Soften { after, floor, recover }`、`Force` | enum 列表；软化**按链**记录、成功提交后回升 |
| 调度 `Schedule` | `RoundRobin`、`Sequential` | enum；取代 `serial` 旗标 |
| 分辨率 | 由 `Space` 决定（连续：全原子；格：重原子 + 装饰） | 装饰是 `Space::finish_chain` |

三个入口成为预设，公开 API 保持（除下面列出的删除）：

- `CbmcGrow(prior)` ≡ `Grow<Continuum, OverlapField>` + `Rosenbluth { Sample(12), 2.0 }` +
  `[Retract, Soften, Force]` + `RoundRobin`；
- `LatticeGrow(prior)` ≡ `Grow<Diamond, Occupancy>` + `Rosenbluth { Enumerate, ∞ }` +
  `[Retract, Reseed, Soften{guard off}, Force{zigzag}]` + `RoundRobin`（**改**：原为顺序）；
- **新** `WalkGrow(prior)` ≡ `Grow<Continuum, NoExclusion>` + `Rosenbluth { Sample(1), 0 }` +
  `[]` + `RoundRobin`——Auhl 路线的第一步：理想链，无排除体积，统计完全由先验决定。

可观测行为：`WalkGrow` 在熔体密度下的 C_n 与稀释极限相同（等于先验的 C∞），R_g 分布
与孤立链一致——"生长偏置导致收缩"由构造消失；`LatticeGrow` 获得每链哈希流与轮转调度，
"删一条链不改另一条链"对格相也成立；`CbmcGrow` 的软化局部化：一条被约束楔死的链不再让
全盒子缩核。

## Non-goals

- 不实现 `Selector::Recoil`（l > 1 feeler）与 `Selector::Bridge`（高斯闭合修饰）：enum 只含
  已实现的变体，两者分别由后续 spec（recoil / ring-closure）加入。
- 不做格上 MC 修复（精修族 R3）、不做格点近似装饰（`topology-model`）：`decorate.rs`
  原样保留，作为 `Diamond::finish_chain`。
- 不改 `InternalTree`、不改先验类型；不改环 / 锚定的拒绝（`topology-model`）。
- 不做 FCC / 立方等第二种格：`Space` trait 是扩展点，但 v1 只有 `Diamond`（earn-complexity）。
- 删除 `relax`（尾部重长）：默认模式下几乎不触发且判据不自洽（评审 A3 / B2）；其职责
  由 `DgRefine`（L4–L5）与 R3（统计）承担。

## Public surface

**Rust 新增**

- `src/grow/space/mod.rs`：`pub trait Space: Send + Sync { type Site: Copy + Send;
  type Var: Copy + Send;  // 连续：扭转角 F；格：RIS 态 i8
  fn n_steps(&self, species) -> usize;
  fn seed(&self, chain: &ChainRef<Self>, draw: &mut Draw) -> Vec<Candidate<Self>>;
  fn step(&self, chain: &ChainRef<Self>, k: usize, draw: &mut Draw, trials: Trials) -> Vec<Candidate<Self>>;
  fn to_continuum(&self, site: &Self::Site) -> [F; 3];
  fn finish_chain(&self, chain: &ChainRef<Self>, coords: &mut [[F; 3]]); }`；
  `pub struct Candidate<S: Space> { pub atoms: Vec<(usize, S::Site)>, pub var: Option<S::Var>, pub prior_weight: F }`。
  **生长单元由 `Space` 定义**（架构师：这不是纯合并而是新设计）：连续空间的一步 = 一个
  `InternalTree` step（种子 3 原子，随后每步一个自由变量 + 其子原子）；金刚石空间的一步 =
  一个主链格点 + 其 RIS 态，全原子在 `finish_chain` 装饰。`Space` 本身不内置四面体键向量——
  `Diamond` 只是一个实现（原则 4；CG 走 `Continuum` + `AnglePrior::Wlc`）。
- `src/grow/chain.rs`：`pub(crate) struct Chain<S: Space> { sites: Vec<S::Site>, vars: Vec<S::Var>,
  stage, visits, deadends, hard_scale: F, .. }`（原 `moves.rs:33-107` 的 `Chain / Species /
  Proposal / Trial`，对空间泛型）+ 哈希流 `stream()`。
- `src/grow/space/continuum.rs`：`pub struct Continuum`（持 `InternalTree` + 角先验 +
  `Seeding { Uniform | VoidBiased }`；`Site = [F; 3]`）。
  `src/grow/space/diamond/mod.rs`：`pub struct Diamond`（`DiamondLattice` + `Backbone`；
  `Site = [i64; 3]`；**上下文无关**：只收 origin / lengths / pbc 数组，永不收 `&PackContext`，
  与今天 `saw.rs:13-18` 一致）；`src/grow/space/diamond/decorate.rs`（原 `lattice/decorate.rs`
  原样搬入，`finish_chain` = `decorate_chain`）。
- `src/grow/exclusion/mod.rs`：`pub trait ExcludedVolume<S: Space>: Send { fn probe(&self,
  atom: usize, site: &S::Site, excl: &[u32], scale: F, soft_shell: F, cap: F) -> Probe;
  fn insert(&mut self, atom: usize, site: &S::Site); fn remove(&mut self, atom: usize);
  fn nearest(&self, atom: usize, site: &S::Site, excl: &[u32]) -> F;
  fn empty_point(&self, u: [F; 4]) -> Option<S::Site> { None } }`；
  `exclusion/field.rs`（`OverlapField`，从 `field.rs` 搬入，377 行不动）、`exclusion/occupancy.rs`
  （`Occupancy`，原 `SawField`，`saw.rs:94-124`）、`exclusion/none.rs`（`NoExclusion`：`probe` 恒
  `Room(0)`，`nearest` 恒 `+inf`——`Force` 的最坏距离排序在它上面仍良定义）。
- `src/grow/select.rs`：`pub enum Trials { Sample(usize), Enumerate }`；
  `pub struct Rosenbluth { pub trials: Trials, pub beta: F }`（`pub type Selector = Rosenbluth`，
  后续 spec 扩成 enum 时不破坏调用方）。
- `src/grow/escape.rs`：`pub enum Escape { Retract { depth: RetractDepth }, Reseed { void_bias: bool },
  Soften { after: usize, floor: F, recover_after: usize }, Force }`；
  `pub enum RetractDepth { Fixed(usize), Exponential { base: usize, every: usize, cap: u32 } }`；
  软化状态存在 `Chain.hard_scale`（**按链**），成功提交 `recover_after` 步后按 `1/0.97` 回升至 1。
- `src/grow/schedule.rs`：`pub enum Schedule { RoundRobin, Sequential }`。
- `src/grow/driver.rs`（重写，≤ 400 行：轮循环 + 死路梯调度）、`src/grow/commit.rs`
  （提交 + 冲突回退，原 `moves.rs:338-391`）、`src/grow/escape.rs`（回撤 / 重播种 / 软化 /
  强制，原 `retract / force_place`）：`pub struct Grow<S, X> { space: S, exclusion: X,
  prior: TorsionPrior, selector: Rosenbluth, escape: Vec<Escape>, schedule: Schedule, .. }`，
  `impl<S: Space, X: ExcludedVolume<S>> Stage for Grow<S, X>`。轮快照语义与
  `stream(seed, mol, stage, visit, salt)` 契约不变，对两种空间统一生效。
- 预设入口拆到 `src/grow/entry/{cbmc,lattice,walk}.rs`（`CbmcGrow` / `LatticeGrow` / `WalkGrow`，
  各 `PackEngine + StageFactory`）；`space/`、`exclusion/`、`driver.rs`、`commit.rs`、`escape.rs`
  **永不** import `crate::entry`（今天 `lattice/mod.rs:28,312` 把 solver 与入口装在一个文件里，
  本 spec 拆开）。`WalkGrow`：`new(prior)`、`with_angle_prior`。
  排除表在 `Target.special_bonds`（只影响后续精修的排除语义，生长期无排除）。`WalkGrow` 的产物按构造有重叠：
  `StageOutcome.converged == false`、`fdist > 0` 是诚实报告，**不**特判。
- `src/grow/config.rs`：`GrowConfig` 增加 `selector / escape / schedule / initial_hard_scale`，
  删除 `relax_every / relax_window / serial / void_bias`（后者进 `Escape::Reseed`）。
  `with_initial_hard_scale(s)`（默认 1.0；dg-refine 的熔体配方用 0.5）。

**Rust 删除**：`src/grow/moves.rs`（拆入 chain.rs / commit.rs / escape.rs / space/continuum.rs）、
`src/grow/lattice/{mod,saw,decorate}.rs`（拆入 entry/lattice.rs、space/diamond/、
exclusion/occupancy.rs）、`src/grow/field.rs`（搬入 exclusion/field.rs）、`src/grow/entry.rs`
（→ entry/cbmc.rs）、`GrowConfig::with_relax / with_serial / with_void_bias`、
`CbmcGrow::with_relax / with_serial / with_void_bias`、`LatticeConfig::with_max_backtrack /
with_max_reseed`（→ `Escape` 列表）。不留别名。

**Python**：`WalkGrow` pyclass；`CbmcGrow` 删除 `with_relax / with_serial / with_void_bias`，
新增 `with_schedule("round_robin" | "sequential")`、`with_reseed_void_bias(bool)`、
`with_initial_hard_scale(f)`；`LatticeGrow` 新增 `with_schedule`；`.pyi` / `_protocols` 同步。
**CLI**：无。

## Module placement

架构师 review-mode 意见（2026-09-02）已合并；按关注点拆分，不按行数拆分
（`moves.rs` 今天的拆分自述"按预算不按语义"，是本 spec 要纠正的）。

| 文件 | 内容 | 预算 |
|---|---|---|
| `src/grow/mod.rs` | 门面 + `validate_template` / `validate_grow_cell` | ≤ 150 |
| `src/grow/config.rs`（叶子，既有） | `GrowConfig` + 吸收 `Selector / Escape / Schedule / Seeding` 数据 | ≤ 400 |
| `src/grow/prior.rs`（叶子，既有） | 不变 | — |
| `src/grow/select.rs` / `escape.rs` / `schedule.rs`（叶子） | enum 定义；只 import molrs | ≤ 150 各 |
| `src/grow/chain.rs` | `Chain<S>` / `Species` / `Candidate` + 哈希流 | ≤ 200 |
| `src/grow/driver.rs`（重写） | 泛型轮循环 | ≤ 400 |
| `src/grow/commit.rs` | 提交与冲突回退 | ≤ 200 |
| `src/grow/escape_ops.rs` | `retract / reseed / soften / force` 的执行（`escape.rs` 只放数据） | ≤ 300 |
| `src/grow/space/{mod,continuum}.rs`、`space/diamond/{mod,decorate}.rs` | `Space` 与两个实现；装饰随格搬 | ≤ 300 各（decorate 393 不动） |
| `src/grow/exclusion/{mod,field,occupancy,none}.rs` | `ExcludedVolume` 与三个实现 | field 377 不动 |
| `src/grow/entry/{cbmc,lattice,walk}.rs` | 三个预设（`PackEngine + StageFactory`） | ≤ 300 各 |
| `src/grow/internal.rs`（600，既有） | 不动（键图辅助已由 stage-pipeline 提升到 `src/topology.rs`） | — |

- **依赖方向**：`entry/* → driver → space / exclusion / commit / escape_ops → chain → config /
  prior / internal / topology`；`space`、`exclusion`、`driver` 不 import `crate::entry`（验收 grep）。
- **叶子纪律**：`config.rs` 只能 import `prior / select / escape / schedule`；四者不得 import
  `space` / `exclusion` / `driver`（否则 entry ↔ grow 成环，与今天 `config.rs:1-6` 的纪律同）。
- `void_bias` 不留在共享 config：它是 `Continuum` 空间的播种选项（`Seeding`）加 `Escape::Reseed`
  的参数——格空间的播种是随机格点。
- **调度与选择器的相容性**：轮快照契约（提议对轮初快照、串行提交）是驱动的定义；v1 的
  `Rosenbluth` 与它相容。未来探测活场的选择器（recoil feeler、bridge 看另一头）**只能**与
  `Schedule::Sequential` 配对，其它组合 `Grow::validate` 具名拒绝——写进 `select.rs` 文档，
  后续 spec 不得绕过。
- 命名：`Occupancy` 取代 `SawField`（SAW 是走法名不是场名）；`NoExclusion` 不叫 `Ideal`
  （理想是统计结果不是排除体积模式）；`WalkGrow` 匹配 `*Grow` 入口模式；无 "packmol"。
- Python：三个预设 pyclass 1:1；`Schedule` 作为 pyclass enum（`python/src/grow.rs` 的
  `TorsionPrior` 模式）；`Grow<S, X>` / `Space` / `ExcludedVolume` / `Escape` Rust-only，
  `Escape` 通过预设的 builders（`with_retract / with_soften / ...`）间接配置。

## Numerical contract

- **CbmcGrow 守门**：~~对原集成 fixture（`grow_pack` 8 × 12 珠、双物种、KG 熔体、约束算例）
  于重构前记录坐标哈希、重构后逐位相同~~——2026-09-29 撤：这些 fixture 已于 2026-09-20 删除，
  且 conventions 禁止 golden / 逐位连续性测试。改为：`src/grow/tests/` 全绿 + 同种子确定性单测
  （形式同 `gencan/entry.rs::gencan_entry_is_deterministic`）。RNG 契约不变：每步的抽样次数与顺序保持（种子 3 + 3 uniform，
  每 trial 1 torsion + 角先验 draws，选择 1 uniform）。
- **LatticeGrow 统计守门**（顺序 → 轮转、全局 RNG → 哈希流，数值必变）：
  原集成测试 `lattice_grow_bead_chain_constructive`（2026-09-20 删除）的 `fdist == 0`、`softened == 0`、
  键长逐位模板三条性质改由 `src/grow/lattice/` 模块内单测承担；稀释极限 RIS 走法的 C_n 与 `prior_ris_calibrated_c_inf` 同容差；新增
  "删一条链不改另一条链"对格相成立。
- **WalkGrow**：稀释极限与 ρ* = 0.85（KG 珠链）两种密度下 C_n 逐位相同（无排除体积 ⇒
  与密度无关，这是可断言的构造性质）；RIS 先验下 C_n = C∞ ± 0.3；`converged == false`
  且 `fdist > 0`（按构造有重叠，如实报告，不特判）。
- **软化局部化**：双物种算例中一个物种带不可满足约束 → 只有该物种的链 `hard_scale < 1`，
  另一物种的最小分子间距 ≥ tolerance − 1e-9。
- 无 ff；`softened` 语义不变（每次按链缩核计一次，`Force` 计一次）。

## Test plan

- `src/grow/tests/{internal,field,prior,entry,driver}.rs`：类型名与 builder 调整；删除 relax 相关断言。
  ~~`grow_cg_kremer_grest_c_inf` 保持 ±10%~~（该 KG 夹具于 2026-09-20 删除，C∞ 由
  `prior.rs::prior_ris_calibrated_c_inf` 承担）。
- 新增属主模块内单测（`src/grow/tests/` 下新子模块，crate 无 `tests/` 目录）：
  1. ~~CbmcGrow 逐位哈希（四个 fixture）~~——撤：禁止 golden，改为同种子确定性；
  2. LatticeGrow 流独立性 + 统计三断言；
  3. WalkGrow 密度无关 C_n（逐位）与 RIS C∞ 命中；
  4. 软化局部化（双物种 + 不可满足约束）；
  5. `Schedule::Sequential` 在两种空间下都可运行且 `softened` 如实；
  6. `Rosenbluth { Enumerate }` 在连续空间报具名错误（无意义组合，`Grow::validate`）；
  7. `Escape` 列表为空时死路 → `converged == false` 且具名 `StageOutcome`（WalkGrow 永不死路，
     CbmcGrow 空梯用于测试）。
- GENCAN 路径（五个 `--example pack_<name>` 程序）不受影响（纯 grow 路径；原 `examples_batch` 已删）。
- `python/tests/test_grow.py`：`WalkGrow` 冒烟、`with_schedule` 校验、被删 builder 不存在。

## Doc plan

- `docs/python/guide/growth.md`：新节「Six axes of growth」（表 + 三个预设的元组）、
  「WalkGrow: ideal chains for melts」（与 dg-refine 配方衔接）；「Melt density」节改写
  （lattice 不再是熔体唯一路径）。
- `docs/python/api-reference.md`（WalkGrow、builder 增删）、`docs/architecture.md` 与
  CLAUDE.md 架构表（`src/grow/` 行按新目录重写）。
- `src/grow/mod.rs` 模块文档：六轴、预设元组、软化按链语义。

## Risks / open questions

1. **逐位守门可能因抽样顺序改变而失败**：若统一驱动必须改变某一步的 RNG 消费顺序，
   规则是——先证明是顺序而非逻辑差异（同分布检验），再更新哈希并在 spec 落地记录里
   写明原因；不允许静默更新。
2. **泛型单态化的编译时间**：三个预设 × 两 trait 组合，可接受；若 `driver.rs` 因泛型膨胀
   超预算，把提交 / 死路梯抽到非泛型的 `commit.rs`。
3. **格相轮转调度的占据表**：轮转要求每步插入 / 删除，`Occupancy` 已是 O(1) 哈希；
   `forced_zigzag` 在轮转下改为按链 `Force`（只放当前步），语义等价。
4. **`initial_hard_scale < 1` 的 `converged`**：生长阶段的 `converged` 只对声明的 scale
   成立；管线裁决在 scale = 1.0，由 dg-refine 负责把它推到容差——文档写明两个 `converged`
   的含义。
5. **`WalkGrow` 的 1-5 自穿**：无排除体积意味着链可能局部自穿（g±g∓ 类构型）；这是 L4
   缺陷，交给精修（dg-refine 的 `IntraOverlapTerm`）；验收只看统计与成键几何。
6. **对 stage-pipeline 的依赖**：三个预设写在今天的 `PackEngine`（`src/entry/mod.rs:108`）
   上即可落地，不阻塞于 stage-pipeline；接缝升级时预设只改 trait 名与 `stages()` 签名。
   但 `src/topology.rs` 叶子（键图 / 排除表）是两者共用的前置，先落地。
7. **这不是纯重构**（架构师 CRITICAL）：格相从全局 `SmallRng` 改为哈希流、从顺序改为轮转，
   数值必变；连续路径逐位保持、格相统计重钉——两者都在 Numerical contract 里各自钉住，
   落地记录写明格相新的参考哈希。
