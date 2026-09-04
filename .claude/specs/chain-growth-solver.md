---
title: chain-growth-solver — a configurational-bias growth solver ranked alongside gencan
status: approved
created: 2026-08-28
revised: 2026-08-28 — 合并 /mol:litrev 结论与四条设计原则；per-target 方法选择取代全局 with_solver；扭转先验升为承重件
---

# chain-growth-solver

## Summary

molpack 今天只有一个 packing 算法：刚体放置 + GENCAN 下降。本 spec 引入第二个
**平级**的算法——构型偏置链生长——并为此在 `pack_with_report` 中段立一条
`Solver` 接缝。

生长**不是**给刚体 packer 做预处理的组件，也**不是**挂在
`Molpack::with_optimizer` 钩子上的构象搜索（该钩子在熔体密度下已实测无效，见
Domain basis）。它是一个独立求解器：消费同一个 `PackContext`（同样的半径、约束、
cell），被同一个目标函数（`fdist`/`frest`）评判，返回同一个 `PackResult`。
**用哪个算法由用户按 target 声明**（`Target::with_method`），不是全局开关，
更不是 molpack 替用户判断。GENCAN 与生长的关系是**并列**，不是包含。

它解决的是刚体模型解不了的一类问题：高密度聚合物熔体。刚体模型对水、尿素、脂质、
蛋白都是对的——那些分子的形状在 packing 开始前就已确定。对熔体它是错的，而且错在
几何而非数值。

## 设计原则（用户裁定，2026-08-28，约束本 spec 全部内容）

1. **molpack 是纯几何的，不碰力场。** 任何构象先验必须是用户提供的几何数据
   （扭转态权重、目标特征比 C∞、持久长度、模板值），molpack 内部永不由力场推导。
   `ff` 只能是可选增强，永不成为 solver 的依赖。
2. **packer（GENCAN）与 grow 平级。** solver 没有任何理由触碰 pack 内部代码
   （`pgencan` / `run_phase` / `run_iteration`），但两者共享同一套架构与生命周期：
   基础设施段 ①②⑤、`PackContext`、共享目标函数、`PackResult`。
3. **用户自己选择什么 target 用什么方法。** 不替用户判断，不静默退化
   （小分子不自动降为刚体）。不支持的组合报具名错误，不猜。
4. **算法须同时适配 all-atom 与 CG。** 排除深度、角度处理、可旋转键感知等决策
   不得硬编码 all-atom 假设。

## Domain basis

### 一、实测：刚体路径在熔体密度下不收敛

体系：PEO 单链模板由 `AmberPolymerBuilder`/GAFF2 建出，25 条拷贝，
`tolerance = 2.0`，seed 42，单线程。`fdist` 的语义是
`max(rsum_ini² − d²)`（`src/objective.rs:293-296`），单位 Å²；`rsum = 2.0` 时
`fdist = 4.0` 意味着两原子完全重合，收敛需要 `fdist < precision = 0.01`，
即最近对 ≥ 1.9975 Å。

dp = 200（1402 原子/链，35050 原子）的密度天花板：

| ρ (g/cm³) | L (Å) | 结果 |
|---|---|---|
| 0.10 | 154.09 | 收敛，3.8 s（3 loops） |
| 0.20 | 122.30 | 收敛，12.0 s（7 loops） |
| 0.30 | 106.84 | 收敛，27.0 s（14 loops） |
| 0.35 | 101.49 | 收敛，24.2 s（11 loops） |
| 0.40 | 97.07 | 60 loops 未收敛，fdist 停在 2.58 |
| 0.50 | 90.11 | 60 loops 未收敛，fdist 停在 2.97 |

天花板落在 ρ ≈ 0.35，约为 PEO 实际熔体密度的三分之一。生产 workflow 用的正是
ρ = 0.5，报错 `molpack did not converge`。

### 二、实测：in-loop optimizer 在熔体密度下的接受率是 0

dp = 100（702 原子/链，17550 原子），ρ = 1.0（L = 56.8 Å），`max_loops = 15`：

| 模式 | 耗时 | converged | fdist (Å²) | 最小分子间距 | 链 Rg (min/mean/max) |
|---|---|---|---|---|---|
| rigid 对照 | 45.9 s | false | **3.3252** | 0.82 Å | 23.1 / 23.1 / 23.1 |
| `TorsionMc` 25 步 30° | 49.0 s | false | **3.3252** | 0.82 Å | 完全相同 |
| ＋`with_self_avoidance(1.0)` | 52.2 s | false | **3.3252** | 0.82 Å | 完全相同 |
| ＋`with_environment(8.0)` | 103.8 s | false | **3.3252** | 0.82 Å | 完全相同 |
| 1 步 2°＋sa＋env | 68.3 s | false | 3.2703 | 0.85 Å | 23.1371 / 23.1404 / 23.1557 |
| 5 步 5°＋sa＋env | 73.6 s | false | 3.6892 | 0.56 Å | 23.1368 / 23.1441 / 23.1858 |

前四行的 fdist 轨迹**逐位相同**（3.6878 → 3.9288 → 3.7087 → 3.9030 → …），
每条链 Rg 亦逐位相同。感知到 1402 根可旋转键，但没有任何一个构象提议存活：
全部被 `src/optimizer/mod.rs:301` 的 non-harm 门回滚。`with_environment(8.0)`
只让墙钟时间 +126%，结果一字不差。

低密度对照（ρ = 0.3，同一体系）：rigid 与 TorsionMc 都在 3.5 s 收敛，
Rg 同为 23.1382/23.1382/23.1382。**即：packing 容易时不需要改构象，
packing 难时改不动。**

### 三、为什么这个钩子救不回来（四条各自独立的缺陷）

1. **步子太大。** `rotate_around_bond`（`src/optimizer/torsion_mc.rs:227`）旋转
   整个 `bond.downstream`；702 原子链的中段扭转一次搬动约 350 个原子。熔体里任何
   这种 pivot 都必然制造新的重叠。
2. **接受判据是全局、贪心、零温。** `src/optimizer/mod.rs:271/300/301` 对**全系统**
   目标函数取 `f_before`/`f_after`，`f_after > f_before` 即回滚。没有 Metropolis、
   没有退火、没有局域性——互相穿插全程是上坡，零温贪心翻不过去。代价是每 copy
   每轮两次全系统求值，只为拒绝。
3. **MC 自己的能量是瞎的且是 O(N²)。** `self_avoidance_penalty`
   （`torsion_mc.rs:206`）是裸双循环，**无 cell list、无最小镜像**；而
   `self_avoidance_radius` 默认 `0.0`（`torsion_mc.rs:73`）使 `energy()` 直接
   返回 0（`:200`），内部 Metropolis 接受一切提议，等于把随机化后的构象交给外层门
   去拒绝。唯一能看见邻居的 `with_environment` 其筛选是每 copy
   O(n_movable × N_total)（`mod.rs:233`）。
4. **跑的时机不对。** `run_optimizer_bindings` 只在 `is_all` 时调用
   （`src/packer.rs:1297`）。单组分体系只有两个 phase，phase 0 先吃光整个
   `max_loops`。

### 四、为什么这不是调参问题

用 PEO 的 Flory 未扰动比 ⟨R²⟩₀/M = 0.805 Å²·mol/g：

| dp | M (g/mol) | 模板 Rg | 熔体理想 Rg | L@ρ=0.5 | L@ρ=1.0 | ρ=1.0 时每点被几条链的 pervaded volume 覆盖 |
|---|---|---|---|---|---|---|
| 25 | 1103 | 6.3 | 12.2 | 45.1 | 35.8 | 4.1 |
| 50 | 2205 | 11.8 | 17.2 | 56.8 | 45.1 | 5.8 |
| 100 | 4407 | 23.1 | 24.3 | 71.5 | 56.8 | 8.2 |
| 200 | 8829 | **46.0** | **34.4** | 90.2 | 71.6 | **11.6** |

两个结论：

- **模板构象本身不对，且错的方向随 dp 翻转**：dp = 25 过度塌缩（6.3 vs 12.2），
  dp = 200 过度伸展（46.0 vs 34.4）。packer 不能信任调用方交进来的形状。
- **dp = 200 在 ρ = 1.0 下，盒子里每一点被约 11.6 条链的 pervaded volume 覆盖。**
  这就是熔体的定义。要求所有原子对不重叠、而搜索手段是刚体平移旋转的算法，必须把
  约 12 条互相贯穿的线团穿过彼此——这是全局拓扑重排，任何从随机刚体投放出发的
  下坡搜索都到不了。熔体必须被**长**成互相贯穿的，不能被**松弛**过去。

### 五、文献调研（/mol:litrev，2026-08-28，两份并行报告的合并结论）

**5.1 扭转先验是承重件，不是可选项。** 固定四面体键角、扭转独立均匀采样是
自由旋转链，C∞ = (1−cosθ)/(1+cosθ) = **2.00**（精确值）。PEO 目标
C∞ = 0.805 × 44.05 / 6.431 = **5.51**（⟨R²⟩₀/M = 0.805 Å²·mol/g [7,8]）。
均匀采样的 Rg 偏差 = √(2.00/5.51) − 1 = **−40%**。熔体排除体积屏蔽修的是
标度指数不是前置因子：亚链修正 ≈ 0.41/√s，s = 600 时对 R² 只有 ~2% [9]——
2% 的效应补不了 2.75 倍的 C∞ 差距。纯几何解法：**RIS 三态先验**
trans(0°)/gauche±(±120°)，由目标 C∞ 单标量解出 trans 分数：
C∞ = C_FRC·(1+⟨cosφ⟩)/(1−⟨cosφ⟩)，PEO ⇒ ⟨cosφ⟩ = 0.467 ⇒ **p_t = 0.645**。
参考实现先例：polyply 收持久长度、Amorphous Cell 收 RIS 表——**两者的生长步
都不用力场** [1,6]，与原则 1 一致。

**5.2 recoil growth 是 CBMC 的严格推广。** feeler 长度 l = 1 时 recoil growth
就是 CBMC；l = 3–8 在高密度下比 CBMC 快 6–50 倍 [3,4]。回撤上界、fallback 位点
复用未试方向、⟨k⟩（可取分数）调到饥饿域边缘，都有现成配方 [3]。
v1 以 l = 1 落地，但 feeler 深度从第一天就是参数。

**5.3 贪心 Rosenbluth 的偏差要声明，不必修。** 不带链级权重簿记的加权选择产出
P ∝ e^(−βU)/W 的有偏系综，偏差随 N 指数增长；硬核极限下穿过拥挤区的链被系统性
高估 [3]。对**构造**（而非平衡采样）这是可接受的固定畸变——本 solver 交付几何，
不交付系综（见 Design §9），产物统计由验收判据直接测量。PERM 类群体控制与
固定链数、固定盒子不相容 [5]，不采用。

**5.4 packer 的交付边界（Auhl 配方 [10]）。** 熔体制备的经典流程是：按目标
C∞ 生成链 → pre-packing 压低密度涨落 → slow push-off → MD 平衡。packer 的职责
是前两步的产物性质：无接触违反、链统计符合先验、密度均匀
（E(d) = ⟨n²⟩−⟨n⟩² 涨落指标 [10]）；平衡系综是下游 MD 的职责。Auhl 的
push-off 停在 0.8σ——`min_hard_scale` 默认 0.8 与之对齐。

**5.5 CG：角度是软坐标。** Martini 类 CG 的角力常数就是为复现持久长度拟合的，
角度**不能**照抄模板；无二面角的 CG 链，持久长度完全由角度先验控制
（WLC 弯曲先验，gensaw `-bendWLC` 先例 [11]）。排除深度是 per-template 属性：
AA 惯例 1-4（深度 3），CG 惯例 1-2 或 1-3（gensaw 显式暴露 `-12/-123/-1234`）。
CG 验证锚点：Kremer–Grest 熔体 c∞ = 1.76（预测）/ ≈1.7（实测）[10,14]。

**5.6 参考实现的 API 形态。** 方法选择在成熟工具里都是 per-molecule 的
（polyply 的 `[molecule]` 块），不是全局开关；先验都是几何数据
（持久长度 / RIS 表 / 扭转直方图）；重构后处理永远归用户 [1,6]。

**5.7 Theodorou–Suter 原文核验（用户提供 PDF，2026-08-28）。** 二手文献转述的
"结合 Meirovitch scanning 的 lookahead"说法**不成立**：1985 原文的生成方案
（其 Eq. 3）是**单步、单向**的——逐键把 RIS 条件概率乘上长程非键能量增量
（Flory 惯例：1-5 起计）的 Boltzmann 因子后重归一化，原文自述"无法预见尚未
生成链段的重叠"。这正是本 spec Design §4c 的形态：**几何先验 × 排除体积权重的
逐步选择**——T–S 是它的特例（每 RIS 态一个候选、无硬核拒绝、无回撤）。他们
不需要死路逃逸，因为不硬拒绝：初猜能量仍高达 10⁶–10⁸ kcal/mol，重叠靠三段
最小化清除（半半径 soft-sphere → 全半径 → 全 LJ，"blowing up the atomic
radii"——radscale 调度与 Auhl push-off 的先声）。本 spec 用硬核拒绝 + 回撤
换取构造保证；软化兜底后的 grow→gencan 串联即其 stage-1/2 松弛的对应物。
recoil growth（§5.2）仍是唯一有文献支撑的 lookahead 升级通道。
定量锚点：长程偏置生长使链相对无扰尺寸**收缩**——expansion factor
⟨r²⟩^½ = 0.85 ± 0.07、⟨s²⟩^½ = 0.92 ± 0.06（Rg 低约 8%）；松弛期间约 36% 的
骨架键翻转扭转态、回向 RIS 分布（bulk t/g = 1.93 vs RIS 1.84）。
**偏差方向可预期：拥挤偏置使生长链轻度塌缩**，ac-006 的 ±10% 容差与此一致。

参考文献：
[1] Grünewald et al., *Nat. Commun.* 13, 68 (2022), 10.1038/s41467-021-27627-4（polyply）
[2] Siepmann & Frenkel, *Mol. Phys.* 75, 59 (1992), 10.1080/00268979200100061（CBMC）
[3] Consta, Wilding, Frenkel, Alexandrowicz, *J. Chem. Phys.* 110, 3220 (1999)（recoil growth）
[4] Consta, Vlugt, et al., *Mol. Phys.* 97, 1243 (1999), 10.1080/00268979909482926
[5] Grassberger, *Phys. Rev. E* 56, 3682 (1997)（PERM）
[6] Theodorou & Suter, *Macromolecules* 18, 1467 (1985), 10.1021/ma00149a018（原文已读，2026-08-28，见 §5.7）
[7] Fetters et al., *Macromolecules* 27, 4639 (1994), 10.1021/ma00095a001（数值经 [8] 转证，未读原文）
[8] Everaers et al., *Macromolecules* 53, 1901 (2020), 10.1021/acs.macromol.9b02428（PEO：0.805 Å²·mol/g，ρ = 1.060 g/cm³ @ 353 K）
[9] Wittmer et al., arXiv:1107.4454（熔体理想性偏差 c_s ≈ 0.41）
[10] Auhl, Everaers, Grest, Kremer, Plimpton, *J. Chem. Phys.* 119, 12718 (2003), 10.1063/1.1628670
[11] Weismantel et al., *Comput. Phys. Commun.* 270, 108176 (2022)（gensaw）
[12] Perez et al., *J. Chem. Phys.* 128, 234904 (2008)（并发生长先例）
[13] Anderson, Irrgang, Glotzer, arXiv:1509.04692（checkerboard 并行 MC + 哈希 RNG）
[14] Kremer & Grest, *J. Chem. Phys.* 92, 5057 (1990)

## Design

### 1. 接缝：`Solver` 是生命周期契约，方法选择在 `Target` 上

`pack_with_report` 今天是五段，其中 ①②⑤ 是基础设施，③④ 才是"算法"：

| 段 | 内容 | 性质 |
|---|---|---|
| ① | 校验 targets、解析 box/cell（含 density → L）、broadcast 全局约束 | 基础设施 |
| ② | 建 `PackContext`：半径、`Constraints`、`SimBox` + `CellGrid` | 基础设施 |
| ③ | 初始状态 | **算法** |
| ④ | 迭代驱动 | **算法** |
| ⑤ | 重建 `xcart` → `assemble_frame` → `PackResult` | 基础设施 |

新增 `src/solver.rs`：

```rust
/// 一个 packing 求解器。实现者把 PackContext 驱动到可行解，
/// 由同一个目标函数评判。solver 永不调用 pack 内部代码
/// （pgencan / run_phase / run_iteration）——共享的是 ①②⑤ 与评判，
/// 不是彼此的迭代器。
pub trait Solver: Send {
    fn name(&self) -> &'static str;

    /// `sys` 已由 ② 建好。`targets` 是本 solver 负责的 target 子集——
    /// 化学（键图、模板坐标）从这里来，与用户传给 pack() 的是同一份，
    /// 不存在第二份真相。实现者写回 `sys.coor` 与 `x`。
    fn solve(
        &mut self,
        sys: &mut PackContext,
        targets: &[Target],
        x: PlacementsMut<'_>,
        budget: &Budget,
        handlers: &mut [Box<dyn Handler>],
    ) -> SolveOutcome;
}

/// gencan 平坐标向量的类型化视图：per-copy 的 com / euler 访问器。
/// 平向量布局是 gencan 的私产，不泄漏给 Solver 实现者。
pub struct PlacementsMut<'a> { /* over &'a mut [F] */ }

#[non_exhaustive]
pub struct Budget { pub max_loops: usize, pub precision: F }

/// 求解结果。fdist / frest 必须来自收尾时对共享 objective 的一次评估
///（`Constraints` 入口），不得由 solver 自报——OverlapField 与 objective
/// 是两套几何代码，验收用同一把尺。
#[non_exhaustive]
pub struct SolveOutcome { pub converged: bool, pub fdist: F, pub frest: F, pub softened: usize }
```

`softened` 的公开载体：`PackResult` **新增** `pub softened: usize` 字段
（加法扩展；既有四字段不变；gencan 路径恒为 0）——ac-004 的 `softened == 0`
断言由它承载，`StepInfo.radscale`（§4g）只是过程可见性。

**方法选择是 per-target 的**（原则 3；polyply 先例 §5.6）：

```rust
pub enum PackMethod { Gencan, Grow(GrowConfig) }   // 默认 Gencan
impl Target { pub fn with_method(self, m: PackMethod) -> Self }
```

（命名取 `PackMethod` 而非 `Method`：后者进 prelude glob 太易撞名，
`Placement`/`CenteringMode` 是既有先例。）`GrowConfig` 与 `GrowError` 放在
**叶子文件 `src/grow/config.rs`**（不 import target/packer/context），由
grow/mod.rs re-export——否则 target.rs ↔ grow/mod.rs 成环（mod.rs 里的
solver 依赖 `Target`）。`PlacementsMut` 与 `init_xcart_from_x`/`SwapState`
是同一布局的两处编码，短期可接受（ac-008 的写回断言守护），`GencanSolver`
抽取时收敛为一处。

早先草案的全局 `Molpack::with_solver(Box<dyn Solver>)` **取消**：全局开关与
per-target 选择冲突，且诱导"solver 替用户判断"的形态。`Solver` trait 仍是公开
接缝（后续 `GencanSolver` 抽取、外部自定义 solver 都在其上）。

**dispatch（分步方案 B 不变）**：`pack_with_report` 在 ③ 之前按 target 的
`PackMethod` 分组。全 `Gencan` → 今天的内联路径**逐字节不变**（`examples_batch`
守门）；含 `Grow` → 生长阶段先跑（见 §6 混合体系）。把 Packmol 路径抽成
`GencanSolver` 是后续独立 PR（Out of scope）。

**接缝的实际位置（架构预检核实）**：dispatch 点在 `init_frame_constants`
（packer.rs:762）与 `let mut x`（:764）之间。两处前置搬移，均为行为保持的
代码搬移、`examples_batch` 守门：(i) cell 解析（packer.rs:781-797）提到接缝
之前；(ii) `initial()` 内部的 SimBox + CellGrid 安装段（initial.rs:524-581）
抽为 ② 级 helper，`initial()` 与 grow 分支共用——grow 不调用 `initial()`
但 `SolveOutcome` 的共享 objective 评估需要已就位的 grid。①处的全局 RNG
（packer.rs:544）首次消费发生在 `initial()` 内，**grow 分支永不从它抽取**
（grow 用 §2 的哈希流），全-Gencan 路径的随机序列因此不变。

**不退化**（原则 3）：`PackMethod::Grow` 配给原子数 < 3、无键模板、或
`template == None` 的 target → 具名错误（提示改用 `PackMethod::Gencan`），
不静默降为刚体。新增错误以 `PackError` 新 variant 表达（repo 惯例：pre-1.0
直接加 variant，不引入 `#[non_exhaustive]`）。

### 2. 状态契约与 RNG 契约

`PackContext` 已是每 copy 一份构象（`sys.coor`，commit `e9a081f`）。生长
**不需要**新的状态模型：

- `sys.coor[copy]` ← 长出来的构象，按质心重心化
- `x[copy].com` ← 该质心；`x[copy].euler` ← `(0, 0, 0)`（精确，非近似）

后果：**`grow` → `gencan` 串联零转换**——两个平级 solver 可组合。

**RNG 契约（并行预埋，v1 串行即生效）**：随机流按 `(seed, copy, step)` 哈希
成 counter-based 独立流（[13] 的纪律），**不用单一全局流**——一条链的采样消耗
不影响另一条链的流。事后从全局流迁移本身就改数值，所以这条必须从第一天生效。

**轮快照语义（并行预埋，v1 串行即按此定义）**：同一轮（同一步索引）内，所有链的
候选提议与打分针对**轮初快照**的场；提交按轮内固定次序串行进行，提交时只对
"本轮已新提交的原子"做增量硬核复查，冲突者重选。这个语义使未来 rayon 化
（并行提议 + 串行提交）与串行**逐位一致**，是免费的决定，窗口只在现在。
`OverlapField::probe` 保持 `&self`（现状即是），插入/删除才是 `&mut`。

### 3. 约束体系复用，且是硬拒绝

`AtomRestraint::f(&self, x: &[F; 3], scale: F, scale2: F) -> F`
（`src/restraint/mod.rs:47`）**本来就是逐点的**。生长时每个候选原子位置求
`Σ_r r.f(&p, ..)`，**大于 0 直接拒绝**——与硬核同一待遇，全拒则回撤。
于是 `frest == 0` 与 `fdist == 0` 一样是**构造保证**，不是收敛希望
（软偏置保证不了"全部原子在球内"，且 Å² 的违反量与无量纲拥挤惩罚量纲不合，
一个 β 伺候不了两个）。同一批 restraint 对象，GENCAN 用梯度消化，生长用
逐点判定消化。刷子、狭缝、孔道内的链免费得到，不发明第二套约束。

### 4. 算法：模板提供化学与先验挂点，形状与统计由先验决定

**(a) 内坐标分解**（`src/grow/internal.rs`，已写 517 行，复核保留）。对模板键图
从链末端（图直径端点）BFS，每个原子改写为对三个已放置原子的
`(bond, angle, torsion)`。键长、键角、非自由二面角逐字照抄模板；同一根可旋转键
上的其它子原子保留相对代表原子的模板偏移（局部几何与手性精确保留）。可旋转键由
`molrs::perceive::rotatable` 感知（未定级键先当单键，与 `torsion_mc.rs` 同策略）。
原子按**步**分组：一步 = 一个自由变量 + 其后所有已确定原子。

排除表在 `Target.special_bonds`（AA 默认深度 3 即 1-2/1-3/1-4，
CG 惯例 1 或 2，见 §5.5）。**注**：1-4 距离跨自由扭转时由 φ 决定（丁烷 trans
3.9 Å vs cis 2.9 Å），排除 1-4 的理由是"1-4 归扭转先验管辖，硬核会错杀
gauche/cis"，不是"由模板固定"——`internal.rs:66-69` 的注释按此修正。

**(a′) 先验**（`src/grow/prior.rs`，新）。纯几何数据，用户提供（原则 1）：

```rust
pub enum TorsionPrior {
    Uniform,                       // 负对照专用：C∞ = 2.00，熔体 Rg 低 40%（§5.1）
    Template { kappa: F },         // von Mises 围绕模板值
    States(Vec<(F, F)>),           // RIS 离散态 (角度, 权重)
}
impl TorsionPrior {
    /// 三态 t/g± 先验，trans 分数由目标 C∞ 解出（§5.1 的闭式）。
    pub fn three_state_from_c_inf(c_inf: F, theta: F) -> Self;
}
pub enum AnglePrior {
    Template,                      // AA 默认：键角照抄模板
    States(Vec<(F, F)>),
    Wlc { persistence_length: F }, // CG：持久长度完全由角度先验控制（§5.5）
}
```

`GrowConfig` **必须显式给 `torsion_prior`，无默认**——均匀采样定量错误
（§5.1），要用必须写出来。`AnglePrior::Template` 之外的选择使键角成为该步的
采样自由度（CG 路径）。候选按先验采样、按 Rosenbluth 权重 `w ∝ exp(−βU)` 选择
（log-w 累加防下溢）。v1 每 target 一个先验；per-torsion-type 覆盖见 Out of scope。

**(b) 盒子从第 0 个原子起就是最终体积**，由 `Molpack::with_density(ρ)` 算出
（见 §7 归属）。没有压缩阶段，密度直接命中。

**(c) 同步逐步生长。** 所有链在每个步索引上轮转推进（每轮随机化链序，
RNG 按 §2 契约），局部密度均匀上升、没有链被饿死。每步 `n_trials` 个候选
（⟨k⟩ 可为分数并可随填充率自适应上调 [3]），对已放下的一切打分
（`src/grow/field.rs`：cell list + 最小镜像 + 排除表，已写 303 行，复核保留）：
硬核 `r_i + r_j`（molpack 自己的 tolerance）以内**拒绝**；硬核到软壳之间计费；
1-`exclusion_depth` 排除，其外计分——链不能穿过自己。
因为硬核是拒绝不是惩罚，`fdist = 0` 是**构造保证**。

**(d) 死路回撤 = recoil growth 的 l = 1 特例。** 某步全拒 → 回撤 `retract` 步
（默认 10，polyply 先例）重长；fallback 位点复用未试方向；同一步反复死路时回撤
深度指数加深。**feeler 深度 `l` 从第一天就是 `GrowConfig` 参数**，v1 实现
l = 1（≡CBMC），代码结构为 l > 1 留位（§5.2：高密度下 6–50 倍加速的升级通道）。
回撤要求重叠场 O(1) 删除——`OverlapField` 的 `cell_of` 反查链表即为此。

**(e) 边生长边松弛，带 W 守卫。** 每 `relax_every` 步，每条链回撤 `relax_window`
步、对此刻更密的场重长。删除旧尾巴前先复算其 log W_old，新尾巴仅当
log W_new ≥ log W_old 才接受——无条件替换可能把好尾巴换成差尾巴，守卫使
松弛单调不劣化。

**(f) 兜底。** 连续失败超过 `soften_after` 次后逐步收缩硬核，下限
`min_hard_scale`（默认 0.8，与 Auhl push-off 的 0.8σ 对齐，§5.4）；每次收缩
计入 `softened` 并进报告；`softened == 0` 才算 `converged`。

**(g) handler 映射。** 每轮发一个 `StepInfo`：`loop_idx` = 轮号，
`radscale` = 当前 hard_scale（软化对 `ProgressHandler` 直接可见），
fdist/frest 在硬拒绝成立期间恒为 0（文档写明该语义）；最终数字由共享
objective 复算（§1 的 `SolveOutcome` 契约）。`Budget.max_loops` 在生长语义下
是重长事件（死路逃逸 + 松弛轮）的总预算上限，文档写明。

复杂度 O(N_atoms × n_trials)，配 cell list；对照刚体路径的 O(迭代 × 对数)。

### 5. 不依赖 `ff`

`src/grow/` 只用 `molrs::{perceive::rotatable, system, store, types}`，均未被
`ff` 门控。`src/grow/` 与 `src/solver.rs` 进 **default feature**，Python wheel
无需新 feature。先验是几何数据：用户（或外部工具）可以从力场推导权重，但
molpack 的 API 只收数据，不收力场（原则 1）。这与 `src/optimizer/`
（`#![cfg(feature = "ff")]`）形成对照。

### 6. 混合体系：串行组合，复用 fixed 机制

用户可以在同一次 `pack` 里混用方法（原则 3 的完整含义；典型场景：聚合物电解质
= PEO 链 `Grow` + Li⁺/溶剂 `Gencan`）。组合是**串行的两 context 方案**
（架构预检判定：bolt-on 侧，descope 检查点预期不触发）：

1. 盒子按**全部** targets 的总质量定容（§7）；
2. 生长阶段：`GrowthSolver` 只长 `Grow` targets（场里只有它们）；
3. 刚体阶段：每个长出的拷贝以
   `Target::from_coords(实验室坐标).fixed_at(质心)` + 零取向进入第二个
   context——`compcart` 在恒等旋转下逐位复现实验室坐标，整条既有 fixed 管线
   （free/fixed 拆分、objective 的 fixed-pair 跳过）零改动地消化它；GENCAN
   只对 `Gencan` targets 求解；
4. 汇总：新增一个跨两个 context 的位置收集函数，把结果按 targets 声明序
   排回（`positions_in_target_order` 只覆盖单 context）。

前置条件：② 段（PackContext 构建，今天是 220 行内联）抽成可调用两次的函数
——行为保持的代码搬移，`examples_batch` 守门，随 Task 9 落地。
链在刚体阶段**逐位不动**（fixed 语义，ac-013 可直接断言）。
**检查点（Task 9）**：若实测抽取仍侵入过深，混合体系缩为具名错误 + 独立
后续 spec；per-target API 本身不变。

### 7. 盒子与密度的归属

`with_density(ρ)` 加在 **`Molpack`** 上（`with_periodic_box` / `with_cell`
旁边），① 阶段解析为盒长——盒子是基础设施，不是 solver 的参数；由密度定容
需要全部 targets 的总质量，也只有 ① 拿得到。与显式盒子互斥（同时给报错）；
**无默认值**：含 `Grow` target 而既无盒子也无密度 → 报错，不回退。

**质量来源**：默认由元素符号查 `Element::atomic_mass`。`Target::from_coords`
建的 target 元素是 `"X"`（CG 路径正是如此）——为此新增
`Target::with_mass(amu_per_copy)` 覆盖；`with_density` 遇到无法解析质量且
无覆盖的 target → 具名错误，不猜（原则 3）。

### 8. Reuse decision

| 现有能力 | 决定 | 说明 |
|---|---|---|
| `Target` / `Target::resolved_radii` | **reuse** | 物种定义、每原子半径沿用；`with_method` 加在它上面 |
| `AtomRestraint` | **reuse** | 逐点 `f()` 硬拒绝，见 §3 |
| fixed 结构机制（`fixedatom` 等） | **reuse** | 混合体系的刚体阶段，见 §6 |
| `assemble::assemble_frame` | **reuse** | 输出帧装配不改 |
| `PackContext` / `Constraints` | **reuse** | 目标函数、判据不改；`SolveOutcome` 由它复算 |
| `handler::Handler` / `StepInfo` | **reuse** | §4g 的映射 |
| `src/gencan/` | **peer** | 并列，不调用；混合时由 `pack` 编排（§6） |
| `src/optimizer/`（`TorsionMcOptimizer`） | **不动** | 对低密度受约束问题仍有效 |
| `molrs::builder::SelfAvoidingWalk` | **不采用** | 无化学、逐链生长高密度下 DeadEnd |
| `initial.rs` | **peer + 一处抽取** | 生长自带 seeding，不复用其刚体投放；其 SimBox+CellGrid 安装段（:524-581）抽为 ② 级共用 helper（§1 接缝，行为保持） |

### 9. 交付物边界（对齐 §5.3 / §5.4）

本 solver 交付**几何**，不交付**系综**：

- 保证：无接触违反（构造性，molpack tolerance 那把尺）、链统计符合用户先验
  （C∞/Rg/内距曲线可测）、密度均匀（E(d) 涨落指标随报告输出）。
- 不保证：平衡 Boltzmann 系综。贪心 Rosenbluth 的畸变（§5.3）是已声明的、
  由验收判据直接测量的构造属性；平衡是下游 MD 的职责（push-off + 平衡跑）。
  这条边界写进 `src/grow/mod.rs` 模块文档。

### 10. 与在途 spec 的关系

- `pair-loop-context-split`（DRAFT）：生长有自己的 `OverlapField`，不碰 pair
  kernel；`GrowthSolver` 读 `sys.coor` / 写 `x`，先落地的一方决定字段访问写法。
  无冲突，仅顺序耦合。
- `triclinic-cell-downshift`（DRAFT）：`OverlapField` v1 只支持正交盒；三斜落地
  后 `dist2`/`wrap` 改走 `SimBox`（field.rs 内留 TODO 锚点）。

## Files

- `molpack/src/solver.rs` — 新。`Solver` / `PlacementsMut` / `Budget` /
  `SolveOutcome`（后两者 `#[non_exhaustive]`）
- `molpack/src/target.rs` — `PackMethod` 枚举 + `Target::with_method` +
  `Target::with_mass`（§7 质量覆盖）
- `molpack/src/grow/config.rs` — 新，**叶子文件**（不 import
  target/packer/context）：`GrowConfig`（`trials` / `feeler` / `selectivity` /
  `soft_shell` / `retract` / `relax_every` / `relax_window` / `soften_after` /
  `min_hard_scale` / `exclusion_depth` / `torsion_prior`(必填) / `angle_prior`）
  + `GrowError`（internal.rs:24 已引用、此前无家）
- `molpack/src/grow/internal.rs` — **已写（517 行），复核保留**；修正 1-4 注释
  理由（§4a）；`exclusion_depth` 参数化
- `molpack/src/grow/prior.rs` — 新。`TorsionPrior` / `AnglePrior` + C∞ 校准
- `molpack/src/grow/field.rs` — **已写（303 行），复核保留**；`probe` 保持
  `&self`；三斜 TODO 锚点
- `molpack/src/grow/driver.rs` — 新。`GrowthSolver: Solver`：驱动循环、
  seeding、回撤、松弛、软化、写回、收尾 objective 评估（体积预算：mod.rs 若
  兼装 config+driver 预估 800–1100 行超预算，故预拆）
- `molpack/src/grow/mod.rs` — 新，**从零写起**（磁盘上不存在此文件）。模块
  文档（平级理由 + §9 交付边界 + 轮快照语义）+ re-exports，只做门面
- `molpack/src/initial.rs` — 一处抽取：SimBox+CellGrid 安装段（:524-581）
  抽为 ② 级 helper（§1 接缝；行为保持代码搬移）
- `molpack/src/packer.rs` — per-target dispatch（§1，搬移 cell 解析 :781-797
  至接缝前）；`with_density`（§7）；`PackResult` 加 `softened` 字段；混合体系
  两 context 组合 + 跨 context 位置收集（§6）。注：packer.rs 已超文件预算
  （1608 行，既有债务），新增逻辑保持薄壳、重活在 grow/driver.rs
- `molpack/src/lib.rs` — `pub mod grow; pub mod solver;` + prelude re-export
  （`PackMethod` 可进 prelude；`GrowConfig` 等从 `grow::` 取）
- `molpack/python/src/packer.rs` — **类型化** `GrowConfig` / `PackMethod` 绑定
  （不是字符串选择），wheel 不加 feature
- `molpack/tests/grow.rs` — 集成测试（往返、随机 vars 不变量、无重叠、密度、
  先验梯度、CG 合成链、约束、串联、混合、确定性）
- `molpack/examples/pack_peo/` — 评测程序转正（`[[example]]` +
  `required-features = ["io"]`，生长模式不需要 `ff`），量化报告含密度均匀性
  E(d)；`out/` 加入 .gitignore
- `molpack/CLAUDE.md` — 架构表增加 `src/grow/` 与 `src/solver.rs` 两行

## Tasks

1. **Add** `src/solver.rs` seam + `Target::with_method`（`PackMethod`）+
   `pack_with_report` 在 ③ 前 dispatch（762/764 之间）。含两处行为保持搬移：
   cell 解析提前、initial.rs 的 box/grid 安装段抽 helper。全 `Gencan` 走今天
   的内联路径。此步**不新增任何算法**，`examples_batch` 五例必须逐条不变。
   ✅ 2026-08-28（6 tests + fast tier 全绿；examples_batch 通过；fmt/clippy 干净）
2. **Add** `src/grow/internal.rs` 复核。RED 先行——(i) 往返测试（模板扭转值重建
   = 模板坐标，‖Δ‖∞ < 1e-9）；(ii) **随机 vars 不变量测试**：任意自由变量下，
   所有键长、键角、非自由二面角、环闭合键长仍等于模板（往返测试抓不住环键
   被误判为自由变量，这条抓得住）。含分支、含环、无可旋转键三个边界用例。
   修正 1-4 注释；`exclusion_depth` 参数化。
   ✅ 2026-08-28（6 复核测试全绿；1-4 注释已改；`from_frame_with_depth` 落地；
   `template_var` 加断言；`dihedral` 的 −π 边界收口）
3. **Add** `src/grow/field.rs` 复核。RED 先行——插入/删除/再插入幂等、跨周期
   边界最小镜像、`nc < 3` 的 27-cell 去重。
   ✅ 2026-08-28（5 复核测试全绿；`probe` 补上与 `nearest` 一致的自 slot 跳过）
4. **Add** `src/grow/prior.rs`。RED 先行——均匀先验孤立链 C∞ = 2.00 ± 0.1
   （解析精确值，是 NeRF + 采样器的回归基线）；`three_state_from_c_inf`
   校准往返。
   ✅ 2026-08-28（4 测试全绿：校准 p_t=0.645、均匀 C_n=2.00±0.1、全-trans
   即模板、RIS 校准 C_n=5.5±0.3 一次命中）
5. **Add** `src/grow/config.rs` + `src/grow/driver.rs` + `src/grow/mod.rs`
   门面：`GrowthSolver` —— hashed per-(copy, step) RNG、轮快照语义（§2）、
   seeding、log-w Rosenbluth 选择、recoil 参数化回撤（§4d）、W 守卫重长
   （§4e）、软化兜底（§4f）、写回契约（§2）。
   ✅ 2026-08-28（9 driver 测试全绿一次通过：构造性熔体 fdist 严格 0、同 seed
   逐位一致、流独立探测器、双物种；examples_batch 在 run_gencan_stages 抽取后
   仍逐条通过。注：v1 回撤深度指数加深实现为 retract·2^(deadends/4)，失败计数
   仅在成功提交时清零；feeler 参数留待 l>1 实现时加入 GrowConfig。2026-08-29
   按架构复查拆分 driver/moves 两文件；按 numerics 复查修正非周期轴 cell 错位、
   pick_ref 共线守卫、penalty_cap 截断可比性）
6. **Wire** 约束硬拒绝（§3）+ `SolveOutcome` 由共享 objective 复算（§1）。
   ✅ 2026-08-28（RestraintTable 硬拒绝进 score/relax；frest==0 构造性判据过；
   不可满足约束经预算加速软化 + force_place 有限终止；修复 visit 计数活锁——
   提议全拒时也消费 visit，防同流重放）
7. **Wire** handler（§4g）。
   ✅ 2026-08-28（每轮 StepInfo：loop_idx=1-based 轮号、radscale=hard_scale、
   xcart 每轮同步供 XYZHandler；should_stop 生效→converged=false）
8. **Add** `Molpack::with_density`（§7：全 targets 总质量、互斥校验、无默认）。
   ✅ 2026-08-29（4 测试全绿：立方盒解析 1e-9、互斥报错、UnknownMass 具名、
   gencan 同样可用；`Target::with_mass` 落地）
9. **Wire** 混合体系两 context 串行组合（§6）：② 抽为可调用两次的函数、
   长出拷贝以 `from_coords + fixed_at` 进第二 context、跨 context 位置收集。
   **检查点**：抽取实测侵入过深 → 缩为具名错误并另立 spec，其余不受影响。
   ✅ 2026-08-29（检查点未触发：`build_context` 抽取后 examples_batch 通过；
   `src/compose.rs` 落地；测试断言链逐位不动且与纯生长参照逐位一致，一次通过）
10. **Add** `tests/grow.rs` 与 `examples/pack_peo` 转正；AA（PEO）与 CG
    （合成 Kremer–Grest 链，无 io）两份量化报告。
    ✅ 2026-08-29（tests/grow.rs 42 测试 + push-off 不搬移探测器（变异验证）；
    pack_peo 转正含内置 PEO 合成器与 E(d)/内距/链序偏差报告。**注**：KG 熔体
    判据（ac-011）与 AA 熔体判据（ac-004/006(3)/010）的最终断言形式待验收
    重新谈判——见 Testing 附注）
11. **Bind** Python：类型化 `Method` / `GrowConfig` + 先验；wheel 不加 feature。
    ✅ 2026-08-29（typed PackMethod/GrowConfig/TorsionPrior/AnglePrior +
    with_density/with_mass/softened；facade+stubs 同步；python 套件 143 全绿
    含 11 个生长测试，端到端生长可用。环境注记：maturin develop 因 sibling
    molpy 的 molrs git-pin 与 path-pin 冲突而失败（预先存在），wheel 经
    maturin build + uv pip install 验证）
12. **Measure** dp=200×25 / ρ=1.0 耗时与质量四数字（耗时、retractions、
    relax_regrowths、softened）进 commit message。
13. **Document** CLAUDE.md 两行 + `src/grow/mod.rs` 模块文档（平级理由 +
    §9 交付边界 + §5.3 系综声明）。
    ✅ 2026-08-29（CLAUDE.md 架构表 4 行 + 回归命令补 --features io；mod.rs
    交付边界与轮快照声明；docs/python 生长指南页 + API 参考；docs/architecture
    模块图；Python builder docstrings；Rust-only 标记；cargo doc 所触文件零警告）

## Testing

- **内坐标往返** ‖Δ‖∞ < 1e-9 + **随机 vars 不变量**（Task 2；后者是环键误判的
  唯一探测器）。
- **无重叠构造性**：`validation::validate_from_targets` 无违反、
  `PackResult::fdist == 0.0`（严格零）、`softened == 0`、独立复算最小非排除
  间距 ≥ tolerance。
- **密度命中**：相对偏差 < 1e-6（解析换算，防单位错）。
- **链统计阶梯**（先验从弱到强，每级都是上一级的回归门）：
  1. 孤立链、Uniform、无排除：C∞ = 2.00 ± 0.1（解析值）；
  2. 孤立链、Uniform、1-5 硬核：测量并**报告**（文献无先验值，不断言）；
  3. 孤立链、RIS p_t = 0.645：C∞ = 5.5 ± 0.3；
  4. 熔体 dp=200×25、RIS 校准先验：平均 Rg 相对 34.4 Å 偏差 ≤ 10%，且严格
     优于模板的 +34%；各 copy Rg 互不相同（max − min > 1 Å）；
  5. 链序无偏：按生长编号前 5 条与后 5 条的平均 Rg 差 < 5%（最后链被挤压的
     探测器）；
  6. 内距曲线 ⟨R²(s)⟩/s 随报告输出（形状对照 §5.1 的 −0.41/√s，不断言）。
- **CG 合成链**（default feature，无 io/ff）：Kremer–Grest 型珠链模板程序化
  合成，`AnglePrior::Wlc` 校准 → c∞ = 1.76 ± 10%。
- **密度均匀性**：E(d) = ⟨n²⟩ − ⟨n⟩²，d ∈ [2σ, 4σ]，随 `pack_peo` 报告输出
  （对照随机放置基线，不设硬门限）。
- **约束生效**：`InsideSphereRestraint` 算例全部原子在球内且 `frest == 0`
  （构造性，§3）。
- **串联**：生长产物直接喂 GENCAN，零坐标变换、首轮 fdist 不劣化；写回契约
  可断言（euler = (0,0,0)，COM = coor 质心，1e-9）。
- **混合**：PEO(`Grow`) + 小分子(`Gencan`) 算例：链原子在刚体阶段逐位不动，
  整体 validation 无违反（Task 9 检查点若缩容，则断言具名错误）。
- **确定性**：同 seed 两次逐位一致；RNG 流按 (seed, copy, step) 哈希——
  单元测试断言"移除一条链不改变另一条链的提议序列"（轮快照 + 独立流的
  联合探测器，也是并行等价性的地基）。
- **回归**：`examples_batch` 五例在 Task 1 后与全部完成后各跑一次，逐条不变。
- **性能**：dp=200×25、ρ=1.0（35050 原子）单线程完成不 DeadEnd，耗时入
  commit message；对照基线"刚体路径 60 loops 不收敛"。

## Out of scope

- **`GencanSolver` 抽取**。接缝已立；搬迁旧路径是独立 PR。
- **三斜盒**。`OverlapField` v1 正交盒；`triclinic-cell-downshift` 落地后路由
  `SimBox`。
- **`.inp` 脚本关键字**。生长暂只有 Rust / Python API；脚本层走 Packmol 兼容
  扩展位，不发明新配置格式。
- **力场扭转采样**。先验是几何数据（原则 1）；从力场推导权重是外部工具的事。
  若将来做，必须 `ff` 可选，不得让 `src/grow/` 依赖。
- **per-torsion-type 先验覆盖**。v1 每 target 一个先验；共聚物的多类型扭转
  （C∞ 反演非唯一）另立 spec。
- **feeler l > 1 的 recoil 实现**。参数与代码结构 v1 就留位（§4d）；实现与
  调参（⟨k⟩ 饥饿域）等 v1 实测回撤率后再定。
- **并行实现**。rayon 化推迟；但**轮快照语义与 hashed RNG 在 v1 串行版即生效**
  （§2），使未来并行（并行提议 + 串行提交，[13] 的 checkerboard 作参照）逐位
  等价于串行。`OverlapField::probe` 保持 `&self`。
- **修 `TorsionMcOptimizer`**。对低密度受约束问题仍有效；要不要修另开 spec。
- **Martini 排除深度默认值**。文献本轮未确证（1-2 还是 1-3）；`exclusion_depth`
  显式可配，CG 文档给指引不给默认断言。
