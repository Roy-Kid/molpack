---
title: lattice-branch-saw-01-walk — 四面体树在金刚石格上的分叉 SAW
status: code-complete
created: 2026-09-05
grilled: true
chain: lattice-branch-saw (01 of 02; successor 02-docs)
---

# lattice-branch-saw-01-walk — 四面体树在金刚石格上的分叉 SAW

## Summary

`LatticeGrow` 把四面体重原子树（线性链、星、梳；重图度数 `d ∈ 1..=4`）作为金刚石格上的分叉自回避行走一次走完：`k = 1` 仍走同一个 `grow_walk`，每个 `InternalTree` 自由变量只抽一次星形方向，同变量的其余重原子按模板 offset 落在剩余格向上（不第二次抽样）。`Backbone.parent` 是 InternalTree BFS 在 not-H 原子上的投影，不是第二套重图 BFS。装饰挂钩用该 site 已有的 `phi - offset`（不要求 `offset == 0`）；对齐三点 `align = [0,1,2]` 没有挂钩。公开合同是 `LatticeGrow::run` 对四面体星/梳返回 `Ok`，`d > 4` 与孤立重原子仍是 `GrowError::NonTetrahedralTemplate`（文案不再写 “branched staged”）。线性珠链的既有构造性门保持绿色。Python / `docs/python` 留给 `lattice-branch-saw-02-docs`。

## Domain basis

金刚石格（位点整数，单位 `a/4`）。A 点全偶且 \(x+y+z \equiv 0 \pmod{4}\)，键 \(+T_k\)；B 点全奇且 \(\equiv 3 \pmod{4}\)，键 \(-T_k\)：

\[
T = \{(1,1,1),\ (1,-1,-1),\ (-1,1,-1),\ (-1,-1,1)\}.
\]

内联恒等式（硬编码验收）：\(|T_k|^2 = 3\)；\(i \neq j \Rightarrow T_i \cdot T_j = -1 \Rightarrow \theta = \arccos(-1/3)\)。线性延续 trans \(\iff\) 新键等于两步前的键。近邻守卫打开时，非键对至少是第二邻 \(a/\sqrt{2}\)。兄弟子格点彼此是第二邻，守卫不禁止它们。

sp³ 分叉的离散自由度是绕入边 \(P\to B\) 的一次 C3：入向占掉 4 个四面体向中的 1 个，剩下 3 向。重原子度 \(d=2\) 占 1 个（线性 v1）；\(d=3\) 占 2 个；\(d=4\) 占满 3 个；\(d>4\) 无法嵌入。兄弟方向彼此 **不是** 线性链的 t/g±；每个 `InternalTree` 可旋转键 `(p,g)` 只有一个变量，后续孩子带模板 offset。把每个子边独立 Rosenbluth 会让 `SawField` 与 `place_step(vars[v]+offset)` 变成同一取代基的两套构象。

- Auhl et al., *J. Chem. Phys.* **119**, 12718 (2003), arXiv:cond-mat/0306026 — Fig. 7 大尺度烙印；方法适用于 branched/star。\(n=100\) 是 KG 珠轮廓。本段不钉熔体 \(C_n\)。
- Consta et al., *J. Chem. Phys.* **110**, 3220 (1999) — recoil 只写了链；子树弹栈来自占据会计，不把 Consta 当树文献。
- Janse van Rensburg & Rechnitzer, arXiv:1107.2162 — 格多边形 ≠ 开 SAW；环仍 `RingTemplate`。

单位：笛卡尔 Å；格点无量纲 a/4；`track_tweak` rad。

## Design

本段推广既有格相叶子，不新增入口、不新增树类型、不改 `InternalTree` 的代表元选取（`Files` 不含 `internal.rs`）。

### 权威

- **拓扑权威**只在 `Backbone`：`atoms`、`parent: Vec<Option<usize>>`、`children: Vec<Vec<usize>>`、`parent_bond: Vec<F>`（Å）、`mean_bond`（Å）、`align: [usize; 3]`、`hooks`、`follows: Vec<Option<(usize, F)>>`（与 `atoms` 平行，只在 `analyze_backbone` 填充：种子/无 site 为 `None`，否则该重原子 `Site::follows`）。禁止 `WalkTree` / `BranchedBackbone`。`LatticeStage` 不重建 `site_of`。`children[j]` 按主链下标递增。`parent_bond.len() == atoms.len()`，根槽不用。
- **格点权威**只在 `Walk.sites: Vec<[i64; 3]>`，与 `atoms` 同一下标。不变量：根以外，`sites[j]` 是 `sites[parent[j]]` 的四个金刚石邻点之一。禁止从 `wrapped.last()` 或线性 `j-1` 延长。
- `LatticeStage` 把 `&backbone.parent`、`&backbone.children`、`&backbone.follows` 借给 `grow_walk` / `forced_zigzag`。原子数是切片长度。`saw.rs` 不 import `decorate`。

### `analyze_backbone`：InternalTree 投影，不是第二套 BFS

not-H 掩码 **只**在本函数：`element.eq_ignore_ascii_case("h")` 是全原子默认；缺元素列 → `"X"`，CG 走每一个原子。不得拷进 `saw.rs`，不得新增 `Target` 氢 API。

重原子邻接 `adj` **只**用来算度数：

- `d == 0`：具名 `NonTetrahedralTemplate`（孤立 / detached）。禁止从行走里跳过留下洞。
- `d ∈ 1..=4`：合法。
- `d > 4`：具名 `NonTetrahedralTemplate`，文案点名度数；禁止 “branched staged” / “staged for the lattice branch phase”。

根 = InternalTree 顺序（`seed_atoms` 再按 `step_sites` 的 `atom`）里 **第一个 not-H 原子**。`parent[j]` = 该原子沿 InternalTree 父链（site `refs[0]` / 树父）走到的最近重原子在 `atoms` 中的下标；`children` 是该投影的逆。禁止再跑一套重图直径 BFS（两套根不一致时挂钩 refs 检查会静默丢钩或误拒合法树）。少于 4 个重原子的既有拒绝保留。重图必须在该投影下连通；走不到的重原子按孤立处理。环仍由 `topology_for_growth` 报 `RingTemplate`。

`align = [0, 1, 2]`：叶子根的 InternalTree 序保证前三个重原子是长度为 2 的键路径。这三点 **没有挂钩**。

### 挂钩：一个变量一个主链 site，offset 不必为 0

`Site::follows` 的权威在 `InternalTree::from_frame`（本段只读）。可旋转键的第一个 BFS site 是代表元（常为氢，`offset = 0`）；走格的重原子是后续孩子（`offset` 一般非零）。要求 `offset == 0` 会跳过主链扭转或误拒全原子 PEO（v1 已能装线性 AA）。

对每个变量 `v`：

1. `B` = 主链重原子中 `site.follows.0 == v` 的集合。
2. `B` 为空 → `NonTetrahedralTemplate`，点名 offset-0 代表原子。具名，不跳过。
3. 合格 `E ⊂ B`：不在 `align`；`parent` 与 `parent[parent]` 存在；`refs[0] == atoms[parent[j]] && refs[1] == atoms[parent[parent[j]]]`（refs 检查，不是 `j >= 3` 切刀）；`follows = Some((v, offset))`（offset 任意）。
4. `E` 非空：挂钩取 `E` 中最小 BFS 下标，记录 `(j, k, li, v, offset)`，一个 `v` 一条。`vars[v] = wrap_pi(phi - offset)`。
5. `E` 空而 `B` 全在 `align`：不对齐点挂钩。`E` 空且 refs 离主链：具名“参考系离开主链”。

挂钩钉在 `tree.step_var(k) == Some(v)` 的 step。非对齐、`d ≥ 2`、有 site 且 `follows is None`：具名非可旋转内部键。叶子 `d = 1` 且 `follows is None` 合法。

### `grow_walk`：每个变量一次自由决策

签名收 `parent`、`children`、`follows: &[Option<(usize, F)>]`。复用 `DiamondLattice`、`SawField`、`RisWeights`、`g_state`、`T_STEPS`、`Walk`。

**对齐前缀。** 先放 `align` 三点，与今日线性 3-site 起步同构：随机自由 `align[0]`；随机自由格邻 `align[1]`；`align[2]` 一个均匀非回头延续（权重 1，记一帧）。这些方向填进 `align[1]` 的已分配位掩码，**之后**才做挂钩抽取。其余主链原子按主链下标递增、父先子后放置。这不是第二条行走函数。

**前缀之后每个父一张 C3 帧。** 入向 + 已放子女占位。挂钩孩子 = `children[p]` 中 `follows[j] = Some((v, _))` 的最小下标，是剩余四面体向里 **唯一** 的 Rosenbluth 抽取。RIS / 五戊烷用父链键 `(gp→p, p→child)`，不用兄弟键。无祖父则权重 1。同一 `v` 的跟随者、以及非对齐 `follows = None`（`d=1` / 刚性），占剩余邻点中与 `follows.offset` 最近 `wrap_pi` 匹配的那个（用切片上已有 offset，不回模板笛卡尔，不第二次抽）。中心+4 叶星 `n_vars = 0`：没有挂钩孩子；剩下两叶从已放的对齐兄弟按 offset 分配。弹该帧释放它写入的全部格点。禁止 `tried: Vec<Vec<[i64; 3]>>`。`max_backtrack` 每弹一星 +1。

`d=2, k=1` 时该规则退化为今日线性 SAW（同一函数）。`d=3` 占 2/3 向、`d=4` 占 3/3 向仍然成立，但第二、第三槽是确定的。

放置 `j` 永远从 `sites[parent[j]]` 走出，不从 `last()`。

`forced_zigzag` 同样借这三片切片：挂钩孩子在剩余槽里优先 trans（否则 `T_STEPS` 序跳过回头），跟随者按 offset 填；忽略 `free_for` 但仍写入占据；永不失败。禁止线性 `n` 之字。

`LatticeStage` 失败梯不变：守卫 → 无守卫 → 强制嵌入；后两级计 `softened`。

### `decorate_chain`

沿父边重建连续目标：`w[root] = lat.to_continuum(sites[root])`，`w[j]` 从 `w[parent[j]]` 沿 `sites[j]-sites[parent[j]]` 方向按 `parent_bond[j]`（Å）缩放。删除 `windows(2)` / `bonds[j-1]` 路径累加。

相位：`k0 = min(挂钩 step)`（无挂钩则为 `n_steps`）；`place_seed` + `place_step(0..k0)`；用 `coords[atoms[align[i]]]` 与 `w[align[i]]` 求刚体；再解挂钩 step 的 `vars` 后 `place_step`。挂钩 step 的 `step_sites` 不得含 `align` 原子。

目标二面角是 **四个点**：

- 存在 `parent³[j]`：`want = dihedral(w[ggp], w[gp], w[p], w[j])`，`beta` 用 build 系同一四点（第四点 `pos0`）。线性 `j-3..j` 是这条父链的特例；禁止用 BFS 下标 `j-3` 冒充曾祖。
- 否则 trial-at-0：`beta = dihedral(coords[refs[2]], coords[refs[1]], coords[refs[0]], pos0)`，`want = dihedral(coords[refs[2]], coords[refs[1]], coords[refs[0]], wb[j])`。`coords[refs[2]]` 来自 `place_seed` / 前缀 step，**禁止**用 `Walk.sites` 去索引 `refs[2]`（常为氢，没有格点）。

`vars[v] = wrap_pi(phi_lat - offset)`；`track_tweak` 仍只作用在挂钩代表子上。氢与刚性侧基走既有 `InternalTree` 步骤。

格常数仍是全盒一条 `DiamondLattice`（`mean_bond` 按拷贝加权）。禁止按物种各拟合一条格。

### 公开文档（本 crate）

`LatticeGrow` / `lattice/mod.rs` rustdoc 英文合同：tetrahedral heavy tree with degree ≤ 4; linear is the \(d=2, k=1\) degeneracy; degree > 4 is `GrowError::NonTetrahedralTemplate`. 删除 “linear-only / branched staged”。`LatticeConfig::with_track_tweak` 改为 parent-chain tracking（rad）。`GrowError::NonTetrahedralTemplate` rustdoc 同步。`entry.rs` **只改 rustdoc**，不改 `validate_targets` 逻辑。不改 `python/` 与 `docs/python/`（02）。

### Never

- 禁止第二套重图 BFS 当 `parent` 的家。
- 禁止每个子边独立 Rosenbluth。
- 禁止 `WalkTree`、`if linear` 第二条路径、路径 `j-1` 延长。
- 禁止 `offset == 0` 才挂钩。
- 禁止对齐三点挂钩。
- 禁止 `Target` 装饰集 API。
- 禁止把星形 / recoil / parent 抽取只放在入口层测试（现 `src/grow/tests/refusals.rs`）。
- 禁止本段改 `python/` 或 `docs/python/`。

### Reuse decision

- reuse `DiamondLattice`, `SawField`, `RisWeights`, `Walk`, `T_STEPS`, `g_state`
- reuse `InternalTree` / `Site::follows` / `nerf` / `dihedral` / `place_step` — 只读，不改代表元
- reuse `topology_for_growth` / `tree_from_target` / `RingTemplate`
- reuse `GrowError::NonTetrahedralTemplate` — 同一变体，改文案
- reuse `LatticeGrow` / `LatticeStage` / `LatticeConfig`
- reuse `src/grow/tests/mod.rs` `branched_parts`；翻转 `lattice_grow_rejects_branched`
- generalize `Backbone` / `analyze_backbone` / `decorate_chain` — InternalTree 投影 + 按变量挂钩 + 父边重建 `w`
- generalize `grow_walk` / `forced_zigzag` — 借 parent/children/follows；每变量一次决策；子树 recoil
- pattern `src/grow/moves.rs` `retract` — 不调用
- new — `WalkTree` / `Space` / Target 氢 API：不赚

## Files to create or modify

- `src/grow/lattice/decorate.rs`
- `src/grow/lattice/saw.rs`
- `src/grow/lattice/mod.rs`
- `src/grow/lattice/lattice_grow.rs`
- `src/grow/lattice/config.rs`
- `src/grow/config.rs`
- `src/grow/tests/refusals.rs`（当时为集成测试文件；`tests/` 已于 2026-09-20 删除，入口层测试迁入此处）
- ~~`regressions/lattice-branch-saw-01-walk.md` (new)~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）

## Tasks

- [x] Write failing unit tests for `analyze_backbone` / hooks / alignment / degrees (`src/grow/lattice/decorate.rs` `#[cfg(test)] mod tests`)
- [x] Write failing unit tests for `grow_walk` / `forced_zigzag` / one-decision-per-var / k=1 (`src/grow/lattice/saw.rs` `#[cfg(test)] mod tests`)
- [x] Write failing public tests in the entry-level test file (then the integration layer, now `src/grow/tests/refusals.rs`) (replace `lattice_grow_rejects_branched`; keep `lattice_stage_requires_none_guarantees_all` and `lattice_grow_bead_chain_constructive`)
- [x] Generalize `Backbone`, `analyze_backbone`, and `decorate_chain` in `src/grow/lattice/decorate.rs`
- [x] Generalize `grow_walk` and `forced_zigzag` in `src/grow/lattice/saw.rs` (borrowed parent/children/follows; star-decision stack; drop `tried` Vec)
- [x] Wire slices through `LatticeStage` in `src/grow/lattice/mod.rs`; update rustdoc in `entry.rs`, `mod.rs`, `config.rs`, and `GrowError::NonTetrahedralTemplate` in `src/grow/config.rs`
- [x] ~~Add regression example `regressions/lattice-branch-saw-01-walk.md` (public API only; hard-coded goldens, no third-party runtime)~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）
- [x] Verify tetrahedral star/comb complete, `d>4` named reject without "branched staged", AA nonzero-offset hook, linear bead chain still constructive
- [x] Run full check + test suite

## Testing strategy

单元测试在所有者模块内（`conventions.md`）。`src/grow/tests/refusals.rs` 只钉公开 `LatticeGrow::run` / `LatticeStage` 声明。绿条：`cargo test -p molcrafts-molpack --lib -- grow::lattice` 以及 `-- lattice_grow`。

**`decorate.rs` in-module**

- 线性无 element：`parent[j] == Some(j-1)`，`align == [0,1,2]`，hooks 不落在 0/1/2。
- 四面体星（中心 + 4 叶，叶在 `T_STEPS` 方向上、键长 1.53 Å，不是 90° 平面星）：一个节点 `children.len()==4`；`n_vars == 0`（终端键不可旋转）。
- 梳：`branched_parts` 几何，`analyze_backbone` 不报错。
- AA 丁烷：offset-0 site 为 H 时挂钩落在同一 `v` 的重原子上并记录非零 offset。
- `d == 0` / `d > 4` 具名拒绝；文案不含 `branched staged`。
- parent 来自 InternalTree 序的第一个重原子，不另跑直径 BFS。

**`saw.rs` in-module**

- `k = 1`：`n=12`，父边是格邻，相继非回头 \(T_i\cdot T_j = -1\)。
- 星 / 梳：one-draw-per-var 用 **有挂钩孩子** 的夹具（`branched_parts` 梳或带臂的星），不是 `n_vars = 0` 的 5 原子星。刚性 4 叶星只钉占据与 `sites[j]` 邻接父点。
- recoil：弹挂钩决策后跟随者格点不留在 `occ`。
- `forced_zigzag` 永不 `None`，`sites[j]` 邻接 `sites[parent[j]]`。
- 源码无 `tried: Vec<Vec<_>>`。

**入口层（现 `src/grow/tests/refusals.rs`）**

> 2026-09-29 现状：`lattice_grow_rejects_degree_gt_4` 与 `lattice_stage_requires_none_guarantees_all` 在 `src/grow/tests/refusals.rs`；
> star / comb 两条完成测试与 `lattice_grow_bead_chain_constructive` 已随 `tests/` 于 2026-09-20 删除，
> 星 / 梳几何改由属主单测覆盖（`src/grow/lattice/decorate.rs::{analyze_backbone_tetrahedral_star, analyze_backbone_comb_has_branch_children}`、
> `src/grow/lattice/saw.rs::forced_zigzag_embeds_a_star`）。以下为落地时的清单。

- `lattice_grow_tetrahedral_star_completes`：中心 + 4 叶，2 拷贝，盒 20 Å → `Ok`，`natoms() == 10`。
- `lattice_grow_tetrahedral_comb_completes`：`branched_parts`，2 拷贝 → `Ok`，`natoms() == 24`。不要求 `converged`。
- `lattice_grow_rejects_degree_gt_4`：文案不含 `branched staged`。
- 删除 `lattice_grow_rejects_branched`。
- 保留 `lattice_grow_bead_chain_constructive`、`lattice_stage_requires_none_guarantees_all`。

~~**Regression** `regressions/lattice-branch-saw-01-walk.md`：5 原子四面体星，键长字面量 1.53，2 拷贝，seed 7，盒 20 Å，`Ok`，`natoms == 10`，中心–叶 1.53 ± 1e-6。日期 2026-09-05；无第三方运行时。~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）

## Out of scope

- `python/` 与 `docs/python/`（后继 `lattice-branch-saw-02-docs`，含 `test_lattice_refuses_star`）
- `src/grow/internal.rs` 代表元选取
- `Space` / grow-axes 合并 / 格上 MC 修复 / 环闭合 / 第二种格子
- `WalkTree`、Target 氢 API、新 `GrowError` 变体
- 熔体 \(C_n\) / \(P_2\) / Auhl \(n>100\) 统计
- 按物种拟合格常数

---

## 已被取代（2026-09-07，`dad53fe` / `7d69e8a`）

本 spec 写就时，装饰按模板内坐标重建骨架，`track_tweak` 是把重建结果拉回格点的旋钮。之后确立的原则是 **模板只提供拓扑，不提供几何**：骨架重原子直接就是走法选中的格点，不再重建。

因此本文以下条目已不再描述代码：

- §单位 的 `track_tweak` rad，与 §`decorate_chain` 的「`track_tweak` 仍只作用在挂钩代表子上」——
  `LatticeConfig::track_tweak` / `LatticeGrow::with_track_tweak` 已删除（重建不存在了，旋钮没有指涉对象）。
- §公开文档 的「`LatticeConfig::with_track_tweak` 改为 parent-chain tracking（rad）」——同上。
- 骨架键长/键角不再是模板的：键长是格点步长（`DiamondLattice::fit` 按模板平均骨架键长定，
  再被盒子的整除条件推开百分之几），键角是格点的四面体角（仅当三轴取整方式相同时精确）。

仍然成立：分叉 SAW 的走法、`Backbone` 的 InternalTree 投影、按变量挂钩、对齐三点不挂钩、
一盒一格、`GrowError::NonTetrahedralTemplate` 的公开合同。
