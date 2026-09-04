---
slug: dg-refine
criteria:
  - id: ac-001
    summary: DgRefine 是 Stage 接缝上的精修阶段，只依赖数值内核、无 ff、无刚体驱动
    type: code
    pass_when: |
      `src/refine/{config,mod,terms,cartesian,minimize,ladder}.rs` 存在；`config.rs` 只 import
      molrs（叶子）；`DgRefine: PackEngine + StageFactory`；`Term` trait 与七个内置项
      （Overlap / IntraOverlap / Bond / Angle13 / Dihedral14 / Chiral / Restraint）定义在 `refine/terms*`。
      `grep -rn 'gencan::solver\|gencan::phases\|gencan::entry\|run_phase\|run_iteration\|pgencan\|gencan::gencan\b' src/refine/` 无命中；
      `grep -rn 'grow::internal\|use crate::grow' src/refine/` 无命中（键图来自 src/topology.rs）；
      `grep -rn 'cfg(feature = "ff")\|molrs::ff\|molrs::optimize' src/refine/` 无命中；
      `grep -rn 'project_cartesian_gradient' src/refine/` 无命中（梯度直接取 gxcar）。
      `src/gencan/mod.rs` 行数不增；`src/objective/` 已拆分且 `examples_batch` 通过（前置）。
      `cargo build -p molcrafts-molpack`（default）通过；wheel feature 集合未变。
    status: pending

  - id: ac-002
    summary: 每个几何项的解析梯度通过有限差分
    type: runtime
    pass_when: |
      tests/refine.rs 对 OverlapTerm、BondTerm、Angle13Term、Dihedral14Term、ChiralTerm、
      RestraintTerm 及 CartesianObjective 总梯度做中心差分（h = 1e-5 Å）比较，
      相对误差 < 1e-6；用例含跨周期边界的成键对与 fixed 邻居。
    status: pending

  - id: ac-003
    summary: 精修消除重叠且刚体 push-off 做不到的算例它做得到
    type: runtime
    pass_when: |
      8 × 12 珠链 / 26 Å 盒的人工重叠初态经 DgRefine 后 `fdist < precision`；
      两条人工自穿链经 IntraOverlapTerm 精修后，各自 `Target.special_bonds` 表外的同分子最小距离 ≥ tolerance − 1e-9
      （两张 Target 各一张表：深度 1 与深度 3）；
      20 × 24 珠链 / 22 Å 盒（复用 lattice 测试的稠密算例）经 DgRefine 末级
      `fdist ≤ 0.1 × 初态 fdist` 且严格小于同 max_loops 下 `GenCanPack::seeded_from` 的 fdist；
      两者的 `fdist` 都来自管线末尾共享 objective 在 scale = 1.0 的评估。
    status: pending

  - id: ac-004
    summary: 不改拓扑、大尺度统计守恒、成键几何在容差内
    type: scientific
    pass_when: |
      对 `CbmcGrow`（或 grow-axes 的 `WalkGrow`）产物做全阶梯精修：键图逐位不变；
      每链 R_g 相对变化 ≤ 3%；内距曲线 ⟨R²(s)⟩/s 在 s ≥ 50 处相对变化 ≤ 5%；
      精修后 1-2 距离相对模板偏差 ≤ bond_tolerance（默认 0.10），1-3 距离 ≤ 1.5 × 该值。
      在珠链（default feature）与 PEO dp=100 × 25 / ρ = 1.06（examples/pack_peo refine，io）
      两个算例上成立；后者的每级 fdist、末级最近对分布、R_g 变化写进 spec 落地记录。
    status: pending

  - id: ac-005
    summary: 约束、手性、确定性、中途终止
    type: runtime
    pass_when: |
      `InsideSphereRestraint` 下精修后 `frest == 0.0`；给定签名四元组时精修后符号不变；
      同 seed 两次逐位一致（含 with_anneal > 0）；中途 `should_stop` → 后续级不运行、
      `converged == false`、键长仍在容差内。
    status: pending

  - id: ac-006
    summary: objective.rs 抽取是行为保持搬移
    type: runtime
    pass_when: |
      `bin_xcart_into_cells` 与 `accumulate_cartesian_fg` 抽出后
      `cargo test --release --features io --test examples_batch -- --ignored` 五例通过，
      fast tier 全绿；GenCanPack 路径的坐标逐位不变（tests/pipeline.rs 的单阶段等价测试）。
    status: pending

  - id: ac-007
    summary: Python 镜像与文档同步
    type: code
    pass_when: |
      `DgRefine` pyclass 含全部 builders；`python/tests/test_refine.py` 通过；
      docs/python/guide/growth.md 的 push-off 建议改为 DgRefine，新节「Refining with
      distance geometry」存在；api-reference 与 CLAUDE.md 架构表更新；`cargo doc` 零警告。
    status: pending
