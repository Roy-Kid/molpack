---
slug: triclinic-cell-downshift
criteria:
  - id: ac-001
    summary: 五个官方例子仍然收敛且无 restraint 违反
    type: runtime
    pass_when: |
      固定 seed 下 mixture / interface / bilayer / spherical / solvprotein 五例均收敛
      （`cargo run --release --features io --example pack_<name>`；返回的
      `State::fdist <= precision` 且 `State::frest == 0`）。
      （2026-09-29 改写：原判据的 validation 报告与 examples_batch 集成 harness 已删除，
      裁决改读 `State`。）
      **不要求**最终目标函数值与改动前逐位相同：cell 语义本身在修正之列，
      与被替换实现一致不是目标。
    status: done
  - id: ac-002
    summary: src/cell.rs 已删除且无残留引用
    type: code
    pass_when: |
      molpack/src/cell.rs 不存在；`grep -rn "crate::cell" molpack/src` 无命中；
      cargo build 通过。
    status: done
  - id: ac-003
    summary: 六方晶胞打包在真实最小镜像下满足容差
    type: runtime
    pass_when: |
      六方胞（a=b，γ=120°）中打包溶剂到 tolerance t 后，用独立的 27 镜像
      暴力扫描（不复用 packer 自身的 cell list）检查所有分子间原子对，
      最小距离 >= t - 1e-9。强倾斜三斜胞同样通过。
    status: done
  - id: ac-004
    summary: 三斜胞密度达到请求值，且优于正交外接盒裁剪的变通法
    type: scientific
    pass_when: |
      在体积 V 的六方胞中放入 N 个分子，打包成功后 N/V 等于请求密度；
      作为对照，正交外接盒打包 + 裁剪在同一晶格上得到的密度严格更低，
      差值在图中报出（这是论文的 Case 1 图）。
    status: pending
  - id: ac-005
    summary: 混合周期性 slab 不再停滞
    type: runtime
    pass_when: |
      xy 周期、z 由分数坐标 InsideCell 约束的 slab 体系收敛，
      迭代数与墙钟时间同量级于等原子数的全周期体系（不出现现有
      with_periodic_box 强制全轴周期导致的停滞）。
    status: done
  - id: ac-006
    summary: 周期轴上的半空间约束被显式拒绝
    type: code
    pass_when: |
      AbovePlane / BelowPlane 的法向在任一周期晶格方向上分量非零时，
      构建阶段返回 PackError，错误信息点名该轴；
      法向完全落在非周期方向时正常工作。
    status: done
  - id: ac-007
    summary: 固定分子可以跨周期边界
    type: runtime
    pass_when: |
      一个跨越周期面的固定 slab 作为 fixed target 时，打包正常进行，
      且该 slab 与溶剂之间的最小距离在 27 镜像暴力扫描下满足容差
      （Packmol 在此情形直接报错拒绝）。
    status: done
  - id: ac-008
    summary: 小 celldim 下 stencil 既不丢邻居也不双计
    type: code
    pass_when: |
      当某轴的 celldim 为 1 或 2 时（三斜胞与薄 slab 下很常见）：
      neighbor_cells_f/g 预计算表中每个无序 cell 对恰出现一次；
      且非周期轴上 celldim==2 时两个 cell 互相在对方的全壳里
      （被替换实现在此丢邻居，见 molrs cell-grid-api AC-002）。
    status: done
  - id: ac-009
    summary: ~~端到端性能灾难告警~~
    type: performance
    pass_when: |
      STRUCK 2026-09-29: benches/ 已于 2026-09-20 删除且无替代（新度量系统待 spec），以下为历史记录。
      benches/pack_end_to_end 与 benches/pair_kernel 的正交基准
      <= 改动前基线 * 1.10；新增的三斜变体作为长期基线记录首次数值。
      本条只作灾难告警，不构成性能主张。
      实测（改动前 e5ae159 vs 改动后，同机同 molrs，各两次；节点漂移约 2%）：
      compute_f 4.054/4.133 -> 4.287/4.311 us (+5.0%)；
      compute_fg 5.936/6.044 -> 5.907/5.775 us (-2.3%)；
      pack_end_to_end 985.8 -> 1001.0 us (+1.5%)。均在 1.10 门限内。
      注：首次测得 compute_f +16.5%（超门限），根因是 pbc_constants 每次求值克隆
      SimBox、以及丢失了"无周期性直接返回"的短路；两者已修（959f652）。
    status: done
  - id: ac-010
    summary: 质量闸
    type: runtime
    pass_when: |
      cargo fmt --all --check、cargo clippy --all-targets -- -D warnings、
      cargo test（含 default 与 rayon 特性）全部 exit 0。
    status: done
---

# Acceptance — triclinic-cell-downshift

## AC-001 — 五个官方例子仍然收敛且 validation 干净

不做向后兼容：被替换的 cell 语义在非周期轴上是错的（见 molrs 侧
`cell-grid-api` 的 AC-001/AC-002），与它逐位一致不是目标，钉住它反而会把缺陷固化。

这条要守的是**能力**不是**数值**：五个例子仍然收敛、约束仍然满足。数值层面的正确性
由 AC-003 的 27 镜像暴力扫描负责，那是独立 oracle，不是旧代码。论文 Case 0 的
兼容性论据同样按这个口径写——"结果等价"指约束满足度与结构统计，不指逐位复现。

## AC-002 — src/cell.rs 已删除且无残留引用

"下沉"必须真的完成。留着私有正交实现意味着两套语义并存，日后必然分叉。

## AC-003 — 六方晶胞打包在真实最小镜像下满足容差

用**独立**的 27 镜像暴力扫描当 oracle，而不是 packer 自己的 cell list——
否则一旦 cell list 的三斜逻辑出错，验证会和被验证对象一起错。

## AC-004 — 三斜胞密度达到请求值，且优于正交外接盒裁剪的变通法

这是 Case 1 的科学论点：Packmol 只能用正交外接盒，SCM/AMS 文档也承认非正交下
"密度通常低于请求值"。这条把"能/不能"变成一个可以画在图里的数。

## AC-005 — 混合周期性 slab 不再停滞

顺手修掉已知的 `with_periodic_box` 强制全轴周期问题。界面体系（xy 周期、z 受限）
是聚合物/电化学最常见的形态，停滞会直接挡住 Case 1 的目标体系。

## AC-006 — 周期轴上的半空间约束被显式拒绝

这是对 m3g/packmol#121 里维护者自问的 "what does 'above plane' mean in the presence
of PBCs?" 的正面回答：无定义就报错，而不是像 Packmol 那样在"第一个盒子"的参考系里
静默求值——那样约束是否满足取决于用户把原点放在哪。

## AC-007 — 固定分子可以跨周期边界

Packmol 在 `initial.f90` 硬性拒绝这种情形。分数坐标 + 按轴钳位/回绕之后，
跨面的固定 slab 可以被正确索引，这是界面体系的实际需求。

## AC-008 — 小 celldim 下 stencil 既不丢邻居也不双计

molpack 挂起过的"ncells < 3 时半壳计数"问题在此一并了结。molrs 侧实测的结论是：
双计并不发生（`nc > cell` 过滤加去重已经挡住），真正的缺陷是**全壳丢邻居**——
非周期轴上只有两个 cell 时，上面那个 cell 的全壳为空。三斜胞与薄 slab 常落在这个
regime，所以两个方向都要钉。

## AC-009 — 端到端性能灾难告警

明确记录：本 spec 不做任何性能主张。时间只在论文里作为兼容性注记出现。
这里的闸门是"没有灾难性退化"。

## AC-010 — 质量闸

项目标准闸门，含 rayon 特性组合。
