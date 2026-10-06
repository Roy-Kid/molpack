---
slug: collective-com-restraints
criteria:
  - id: ac-001
    summary: Atoms 归约路径与现有行为完全一致
    type: code
    pass_when: |
      SiteReduction::Atoms 下，现有全部 collective restraint 测试逐值通过；
      未绑定 collective restraint 的打包结果与改动前逐字节相同。
    status: pending
  - id: ac-002
    summary: Com 归约的前向与梯度散射正确
    type: code
    pass_when: |
      对一个已知拷贝，R_c 等于按权重归一后的加权中心；
      ∂L/∂r_i = w_i · ∂L/∂R_c 与中心差分一致（相对误差 < 1e-6）；
      ComWeights::Mass 在 target 缺质量时于构建阶段报错。
    status: pending
  - id: ac-003
    summary: 倒格矢枚举完备且去重
    type: code
    pass_when: |
      给定 SimBox 与 q_max，枚举出的 Q 集合等于 {2π (h⁻¹)ᵀ n : n ∈ ℤ³\{0}, |q| <= q_max}
      在 q ↔ −q 对称下的商集；正交与强倾斜三斜胞各验一次；
      非周期轴不贡献倒格方向。
    status: pending
  - id: ac-004
    summary: StructureFactor 梯度与中心差分一致
    type: code
    pass_when: |
      StructureFactor 所属模块的 `#[cfg(test)]` 与
      src/restraint/geometric/tests/gradient.rs（经 objective 的中心差分，同
      `collective_restraint_gradient_matches_finite_difference_through_the_objective` 的形式）中对 StructureFactor 项做中心差分，
      在正交胞与三斜胞上相对误差均 < 1e-6；axis 几何同样通过。
    status: pending
  - id: ac-005
    summary: 均匀熔体：小 q 涨落被压低至少一个数量级且不破坏打包
    type: scientific
    pass_when: |
      同一构象池、同一密度、同一 tolerance 下跑两遍（仅 pairwise 目标 vs 加入
      StructureFactor）。加入后最低若干 q 壳层的 S(q) 相对未加入时下降 >= 10 倍；
      两次运行返回的 `State` 均 `fdist <= precision` 且 `frest == 0`（分子间最小距离满足 tolerance；原 validation 报告已删除）。
    status: pending
  - id: ac-006
    summary: 目标剖面：质心分布匹配给定的层状剖面
    type: scientific
    pass_when: |
      以 tabulated 层状剖面为目标打包后，沿层法向的质心分布与目标剖面的
      W₂ 距离低于 spec 中声明的阈值；同时 pair-distance 约束仍满足。
    status: pending
  - id: ac-007
    summary: 取向分布：cos θ 分布匹配目标
    type: scientific
    pass_when: |
      以给定的 cos θ 目标分布打包后，实际分布与目标的 W₂ 距离低于声明阈值；
      报出由此得到的序参量数值。
    status: pending
  - id: ac-008
    summary: 无力场依赖
    type: code
    pass_when: |
      本 spec 新增的全部代码路径在 default 特性（不含 ff）下编译并通过测试；
      `grep -rn "feature = \"ff\"" molpack/src/restraint/collective` 无命中。
    status: pending
  - id: ac-009
    summary: 质量闸
    type: runtime
    pass_when: |
      cargo fmt --all --check、cargo clippy --all-targets -- -D warnings、
      cargo test（default 与 rayon）全部 exit 0；Python 绑定的对应测试通过。
    status: pending
---

# Acceptance — collective-com-restraints

## AC-001 — Atoms 归约路径与现有行为完全一致

归约层是插在几何前面的新一层，默认路径必须是恒等变换。这条保证新机制不会悄悄改变
已有的 profile 型约束的数值。

## AC-002 — Com 归约的前向与梯度散射正确

集体约束的梯度是跨整个物种耦合的，散射写错不会让程序崩溃，只会让优化悄悄走偏，
所以必须对着中心差分钉死。

## AC-003 — 倒格矢枚举完备且去重

周期条件下只有与晶格公度的波矢是允许的，所以 Q 不是建模选择，而是给定 q_max 下的
完整可容许集合。倾斜胞的 (h⁻¹)ᵀ 是最容易写错的地方，必须单独验。

## AC-004 — StructureFactor 梯度与中心差分一致

解析梯度 ∂S/∂R_c = (2/N)[B_q cos(q·R_c) − A_q sin(q·R_c)] q 的实现验证。
三斜胞单独验，理由同 AC-003。

## AC-005 — 均匀熔体：小 q 涨落被压低至少一个数量级且不破坏打包

这是 Auhl–Kremer prepacking 的复现基线，也是本 spec 的核心功能判据。
**两臂必须用同一个构象池**——链内统计完全相同，唯一变量是放置目标，
否则数字无法归因。同时必须证明集体项没有以牺牲堆积质量为代价（`State` 的 `fdist`/`frest` 干净）。

注意这条是"复现已知方法"，不是新贡献；贡献在 AC-006 / AC-007 的非均匀目标上。

## AC-006 — 目标剖面：质心分布匹配给定的层状剖面

把目标从常数推广到任意剖面，这是 prepacking 做不到的部分。层状剖面来自 SCFT 或实验，
作为 tabulated 目标输入。

## AC-007 — 取向分布：cos θ 分布匹配目标

Packmol 的 `constrain_rotations` 只能给硬性角度边界，给不了分布。这条证明取向也是
一个可被匹配的集体统计量。

## AC-008 — 无力场依赖

项目约束：力场不能是任何主线功能的必要条件。本 spec 的全部内容都是纯几何/统计的，
这条把它钉住，防止实现时顺手引入 `ff` 依赖。

## AC-009 — 质量闸

项目标准闸门，含 Python 绑定。
