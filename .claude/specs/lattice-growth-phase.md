# lattice-growth-phase — 金刚石格相生长

状态:v1 LANDED 2026-09-01(线性主链全管线:守卫式 SAW → InternalTree 装饰 →
共享 objective 如实评估 → `GencanPack::with_restart` push-off 链)。仍属 DRAFT 的
部分:分支模板(named-reject 占位)、格上 MC 修复(`with_repair_sweeps` 未落地)、
统计验收全套(内距曲线 / P₂ / 20 种子矩阵)。原型与基准:
`~/work/molpack-bench/diamond_saw.py`,slurm array 1871323。

**v1 落地记录(2026-09-01)**:
- 入口 `LatticeGrow`(`src/grow/lattice/{mod,config,saw,decorate}.rs`),
  `LatticeSolver: Solver`;Python `LatticeGrow` 镜像。
- **原型 bug 修正**:diamond_saw.py 的 trans 判据 `b2 == b0` 比较的是相邻键
  (`cur−prev` vs `q−cur`),在金刚石格上因子格符号交替**永不为真**——原型
  实际是均匀权重。落地版 trans ⇔ 新键 == 两步前的键;并加了 RIS 五戊烷
  排除(相邻 g±∓ 直接不提议)。
- 装饰:格只定扭转序列;逐步在 NeRF 参考系里数值解 vars(主链二面角对
  站点扭转斜率恒为 1,一次试探即精确),再整体刚体对齐到走格前三点。
  站点占据仅是提议过滤器,verdict 一律来自装饰后共享 objective(§占据表)。
- 失败梯:守卫 SAW → 无守卫 SAW → 强制全反式 zigzag,后两级计入
  `softened`,不收敛如实报告。
- 实测(peo-tg dp25 全原子 PEO ×200,ρ=1.1 g/cm³,35,400 原子):0.09 s,
  softened=0,键长逐位模板,Rg 11.96±3.12 Å vs Flory 12.25 Å;fdist≈4.0
  的残余接触(氢 + 装饰漂移)交给 seeded push-off。连续 CBMC 在同一格
  无界研磨。

## Goal

新增一个格相生长求解器:把盒子映射到金刚石格子(整数坐标,单位 a/4),每条链的重原子主链作为格上自回避行走生成——3 个非回头延伸方向精确对应 trans/gauche± RIS 态(权重沿用 `TorsionPrior` 标定),排斥体积 = 位点占据 + 近邻位点排斥(O(1) 哈希),死角回退为弹栈。生成后在格上做带 RIS 权重的 MC 修复(蛇行 + 局部移动)恢复链统计,再"装饰"回连续空间:格上决定的扭转序列喂给现有 `InternalTree`,以模板真实键长键角重建全部原子(含氢)。可观测行为:melt 及以上密度(位点占据 φ ≤ 0.65,≈ 2.8× 熔体)的 grow 目标在秒级完成构造——连续版在 ρ≥0.35 等效(φ≈0.08)即研磨失败(基准:dp25-n200,20 种子,连续版 ρ0.4 0/40,格相 φ0.65 10/10 × 0.1 s)。

## Domain basis(先例与他人的解法)

1. **直系血统:Mattice 学派 2nnd 格子桥**。Rapold & Mattice 1996(2nnd = 金刚石格隔位相连)在格上做 PE/PEO/PS 熔体 MC;Doruker & Mattice 1997(Macromolecules 30, 5520)逆映射回连续原子模型——补回母金刚石格上的中间原子 + EM/短 MD;验收用扭转分布、g(r)、S(q)、溶度参数。该管线 2014-2021 仍在产出(PEO Mw 8000、PS 立构规整)。本 spec = 母格子(不隔位)版本做成快速 builder。
2. **生长收缩是标准病征,修复路线有定论**。逐键生长的稠密链系统性偏紧(T–S 1985:⟨r²⟩^½ ≈ 0.85 无扰值;我们实测 C_n 5.2→3.8@φ0.23,同源同量级)。文献中**无人**做"按密度调 t/g 生长权重"的预补偿(A1 检索空白);被验证的两条路:(i) Meirovitch 扫描式 look-ahead(Amorphous Cell = TS + scanning);(ii) **生成后在格上做带 RIS 权重的 Metropolis 修复**(Mattice 路线:权重进 Hamiltonian 而非生长先验;非局部移动——蛇行、Pakula CMA 协同回路在 φ=1 都能动,我们 φ0.23 成本极低)。Auhl 2003 Fig. 7 是关键约束:**n>100 的大尺度统计在生成期烙印,push-off 只修 n≲100**——统计必须在装饰前的格相修好。
3. **高密度四面体格子的有序化风险**:KSY 1986(Macromolecules 19, 2560)证明半刚性参数区存在 order-disorder 转变——验收必须查取向序参量 P₂。
4. **装饰/弛豫协议成熟**:先只开键合项最小化(backward EM1)→ LJ σ/ε 缩放爬坡(PolymerModeler,推广 T–S scaled-LJ;早期弛豫用 LJ 勿用 exp-6)→ capped-force / 位置约束短 MD。本库对应物:共享 objective + GENCAN push-off,不新造。
5. **验收指标行业标准**:⟨R²(n)⟩/n 内距全曲线(非单点 R_g)+ S(q) + 扭转分布(Auhl;Zhang/Kremer 层级回映同一判据)。

## 拓扑范围(2026-09-01 补充,owner 指令)

- **支链(树状模板)是一等公民**,不是扩展项。依据:`InternalTree` 本就是键图的树分解
  (根取最长最短路端点,`internal_roundtrip_branched` 测试钉住分支往返精度)。格相对应物:
  分叉 SAW——sp³ 分叉点占用 3 个延伸方向中的 2 个(金刚石格天然四面体分叉),树按
  BFS 步序生长,recoil 按子树弹栈(回退分叉点须连带其已长子树)。验收加一个分支模板
  (梳形 / 短支 PE 型)贯穿全部统计与装饰测试。
- **环拓扑(ring)显式 named-reject**(`GrowError::RingTemplate` 之类):格上闭环是独立的
  算法问题(闭合偏置采样),留给后续 spec;拒绝而非静默降级。
- 不做非四面体主链(芳环、sp² 共轭)的格相映射;同样 named rejection。

## Non-goals

- 不做格上完全平衡(CAMC/端桥是下游 MC 的事);格相 MC 修复只到内距曲线达标。
- 不做 FCC/立方等第二种格子(earn-complexity:金刚石 only,直到出现第二个调用方)。
  命名注意:金刚石格 = 两套互穿 FCC,四面体键角与 t/g± 同构都来自金刚石;入口名 `LatticeGrow`(owner 定名,2026-09-01):方法族描述,具体格子是配置/文档事实;
  文献两名混用(diamond lattice / tetrahedral lattice),无正统性约束。
- 不替换连续生长求解器——两者是独立入口(见 engine-entry-split spec)。
- 装饰后的残余重叠不在格相消灭——共享 objective + push-off 是唯一裁决与修复权威。

## Public surface

**注(2026-09-01)**:owner 指令拆分引擎入口(删除大一统 `Molpack`,每算法独立入口共享
生命周期 trait)——本节的挂载点以 engine-entry-split spec 为准:格相生长的入口是
**`LatticeGrow`**,下述 `LatticeConfig` 成为其构造参数;拆分前的过渡挂载点(per-target
方法枚举)已随拆分删除,不再存在。

**Rust**
- 入口 `LatticeGrow`(`src/grow/lattice/lattice_grow.rs`,`PackEngine` 的实现之一)。
- `LatticeConfig`(新叶子 `src/grow/lattice/config.rs`,只 import `prior` + molrs):`new(torsion_prior: TorsionPrior)`(强制,无默认——同 `GrowConfig` 的先验条款)、`with_repair_sweeps(n: usize)`(格上 MC 修复扫数,默认待标定)、`with_occupancy_guard(bool)`(近邻位点排斥开关,默认 on)。
- `LatticeStage`(`src/grow/lattice/mod.rs`),实现 `Stage`(`src/stage.rs`);`grow::lattice::SawSite` 等类型域内命名(避免与 `region.rs` 的 cell/lattice 概念混淆)。
- 拒绝错误:`GrowError::NonTetrahedralTemplate`(模板键图含非四面体主链节点时)。

**Python**
- 入口类 `LatticeGrow`(`python/src/packing_methods.rs` 的 `PyLatticeGrow`);`.pyi` 与 `_protocols` 同步。

**CLI**:无(.inp 不新增关键字;spec principle:配置面最小)。

## Module placement(architect 意见,已采纳)

- 新 `src/grow/lattice/{mod.rs, config.rs, saw.rs, decorate.rs}`,每文件 200-400 行预算。
- 方法选择走调用方挑入口(`LatticeGrow`,或 `Pipeline::with_stage` 组合),**不是** GrowConfig 旗标(避免 per-target 决策被 any() 折叠——`driver.rs` 的 `serial` 旗标即该失误模式,列为本 spec 的顺带清理项:serial 迁往方法级或文档声明其全局语义)。
- 格几何(需要 box)在 `saw.rs`(可 import `context`);config 保持叶子。
- **占据表只是提议过滤器**:最终 `StageOutcome` 与 `State::fdist`/`frest` 一律来自装饰后的 `OverlapField` + 共享 objective(单一权威;`src/stage.rs` 合同不变)。
- `decorate.rs` 复用 `InternalTree::place_step_with_angles` 写 `sys.coor`(格上 t/g± 序列 → 模板真实内坐标重建全原子);`assemble.rs` 不动。

## Numerical contract

- 确定性:每链哈希流(`stream(seed, mol, stage, visit, salt)` 同款),同种子逐位可复现;格相 MC 修复同样按流键控。
- 统计:稀释极限采样 C_n 与 RIS 理论一致(现有 `prior_ris_calibrated_c_inf` 同容差 ±0.3);修复后熔体 φ 的内距曲线 ⟨R²(n)⟩/n 相对目标(短链熔体外推)偏差 ≤10%(对齐 ac-006 容差),全 n 范围;P₂ 取向序 < 0.05。
- 装饰:重建键长/键角与模板逐位一致(1e-9,InternalTree 既有保证);扭转态序列与格上决定一致(t/g± 分类逐链核对)。
- fdist 语义不变:装饰后全尺度评估,残余接触如实计入,收敛判据不放松。

## Test plan

- 属主模块内单测(原计划的集成测试层 `tests/` 已于 2026-09-20 删除;已落地覆盖:`src/grow/lattice/saw.rs::{t_steps_are_tetrahedral, grow_walk_k1_sites_are_neighbours, grow_walk_recoils_out_of_dead_ends, forced_zigzag_embeds_a_star, walk_skips_blocked_half_space}`、`src/grow/lattice/decorate.rs::analyze_backbone_*` / `decorate_chain_seats_every_backbone_atom_when_rooted_at_hydrogen`、`src/grow/tests/refusals.rs::{lattice_grow_rejects_degree_gt_4, lattice_grow_empty_region_is_named, lattice_stage_requires_none_guarantees_all}`;下列统计项仍属 DRAFT):亚格宇称/近邻规则单元测试;SAW 有效性(自回避 + 键角恒 109.47°);占据记账(复用 `field_empty_cell_bookkeeping` 模式);同种子逐位确定性;稀释 C_n 统计;修复后 φ0.23 内距曲线;装饰几何逐位核对;`NonTetrahedralTemplate` / `RingTemplate` 拒绝;**分支模板全链路**(格上分叉生长 + 子树 recoil + 装饰往返,复用 `branched_parts()` 几何)。
- `python/tests/test_grow.py`:`LatticeGrow` builder 链与 `run` 冒烟(natoms、converged);repr 子串。
- 20 种子成功率矩阵(`~/work/molpack-bench` 基建)作为性能验收:φ ∈ {0.23, 0.35, 0.5} 全成,单格 < 5 s(Rust)。
- GENCAN 路径(五个 `--example pack_<name>` 程序)不受影响(纯新增路径;原 `examples_batch` harness 已删)。

## Doc plan

- `docs/architecture.md` 模块表加 `src/grow/lattice/` 行;`docs/python/guide/growth.md` 新节(何时选格相:高密度四面体主链)+ `api-reference.md` 两个新类型;CLAUDE.md 模块表同步。

## Risks / open questions

1. **统计修复的收敛预算**:格上蛇行+局部移动在 φ0.23 修复 −30% R² 收缩需要多少扫?(先例定性支持、无直接数字;用原型标定后定 `with_repair_sweeps` 默认。)若不足,升级为 CMA 协同回路(φ=1 先例)。
2. **p_t(φ) 自标定生长**(文献空白 = 潜在新颖点):格上 0.1 s 重生成使割线迭代标定几乎免费,但单标量只保 C_n 不保内距曲线形状——作为实验分支,与 MC 修复 A/B 后择优;Auhl Fig.7 是判据。
3. **异质主链漂移**:PEO 的 C-O 1.43 Å vs 格键 1.54 Å,装饰链相对格路径累积漂移,格相排斥保证随链长退化——依赖装饰后共享 objective 兜底;若 fdist 残余系统性超标,引入格距 a 按主链平均键长取值。
4. **盒子公度**:`with_density` 解出的 L 与 4×(a/4) 网格不公度——格距按 a′ = L/round(L/a) 微调(格上键长仅入拓扑不入几何,无害);逐轴处理。
5. **KSY 有序化**:P₂ 验收若在某参数区超标,需要记录相图边界并 named-reject 或降 φ。
6. 顺带清理:`GrowConfig.serial`/`void_bias` 旗标的 any() 折叠语义(architect 🔴)——本 spec 落地时一并迁移或文档化。

## 参考

Doruker & Mattice, Macromolecules 30, 5520 (1997);Rapold & Mattice, Macromolecules 29, 2457 (1996);Auhl et al., JCP 119, 12718 (2003)(Fig. 7 生成期烙印;push-off 配方);Theodorou & Suter, Macromolecules 18, 1467 (1985);Kolinski, Skolnick & Yaris, Macromolecules 19, 2560 (1986)(四面体格有序化);Pakula, Macromolecules 20, 679 (1987)(CMA, φ=1);Wall & Erpenbeck, JCP 30, 634 (1959)(enrichment);Consta et al., JCP 110, 3220 (1999)(recoil growth);Wassenaar et al., JCTC 10, 676 (2014)(backward 弛豫协议);Zhang et al., ACS Macro Lett. 3, 198 (2014) & Soft Matter 15, 289 (2019)(层级回映与验收);Dietschreit et al., JCTC 12, 2388 (2016)(化学真实四面体格,Roulattice)。
