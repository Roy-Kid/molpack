---
slug: grow-axes
criteria:
  - id: ac-001
    summary: 一个泛型驱动，六轴各自成型，旧循环与旧旗标消失
    type: code
    pass_when: |
      `src/grow/driver.rs` 定义 `pub struct Grow<S: Space, X: ExcludedVolume<S>>` 并
      `impl Stage`（或 stage-pipeline 未落地时 `impl Solver`）；`src/grow/chain.rs`、`commit.rs`、
      `escape_ops.rs`、`src/grow/space/{mod,continuum}.rs`、`space/diamond/{mod,decorate}.rs`、
      `src/grow/exclusion/{mod,field,occupancy,none}.rs`、`select.rs`、`escape.rs`、`schedule.rs`、
      `src/grow/entry/{cbmc,lattice,walk}.rs` 存在；`src/grow/moves.rs`、`src/grow/field.rs`、
      `src/grow/lattice/` 不存在。
      `grep -rn 'with_relax\|with_serial\|with_void_bias\|SawField\|fn relax\b' src/ python/ docs/` 无命中。
      `grep -rn 'use crate::grow::\(space\|exclusion\|driver\)' src/grow/config.rs src/grow/select.rs src/grow/escape.rs src/grow/schedule.rs` 无命中（叶子纪律）；
      `grep -rn 'use crate::entry' src/grow/space/ src/grow/exclusion/ src/grow/driver.rs src/grow/commit.rs src/grow/escape_ops.rs` 无命中；
      `grep -rn 'PackContext' src/grow/space/diamond/` 无命中（格几何上下文无关）。
      新文件 ≤ 400 行（`exclusion/field.rs` 377、`space/diamond/decorate.rs` 393 原样搬入）；
      `grep -rn 'molrs::ff\|molrs::optimize' src/grow/` 无命中（`ff` 已是纯透传 `molrs/ff`，crate 内无 `cfg(feature = "ff")`，原 cfg 判据改为依赖判据）。
    status: pending

  - id: ac-002
    summary: CbmcGrow 逐位不变
    type: runtime
    pass_when: |
      （2026-09-29 改写：原判据的四个集成 fixture 坐标哈希已随集成层于 2026-09-20 删除，
      且 conventions 禁止 golden / 逐位连续性测试。）
      `src/grow/tests/` 全部单测（`internal` / `field` / `prior` / `entry` / `driver`）在重构后全绿；
      新增一条 CbmcGrow 同种子确定性单测（同 seed 两次运行位置与裁决逐位相同，
      形式同 `gencan/gencan_pack.rs::gencan_pack_is_deterministic`）通过。
    status: pending

  - id: ac-003
    summary: LatticeGrow 获得哈希流与轮转调度，统计断言不变
    type: runtime
    pass_when: |
      `lattice_grow_bead_chain_constructive` 的 `fdist == 0.0`、`softened == 0`、
      键长逐位模板三条断言通过；同 seed 逐位一致；新增测试断言删除一条链后其余链
      的坐标逐位不变（流独立）；稀释极限 RIS 走法 C_n = 5.5 ± 0.3。
    status: pending

  - id: ac-004
    summary: WalkGrow 的统计与密度无关且等于先验
    type: scientific
    pass_when: |
      不实施。晶格驱动就是自回避随机行走，不另设理想链入口。
    status: dropped

  - id: ac-005
    summary: 软化按链局部化且可恢复
    type: runtime
    pass_when: |
      双物种算例中物种 A 带不可满足的 InsideSphere 约束、物种 B 无约束：运行结束时
      物种 B 的最小分子间距 ≥ tolerance − 1e-9（从未被软化），`softened > 0` 全部归于 A；
      单物种可行算例中一条链短暂软化后（人为注入）在 `recover_after` 次成功提交后
      `hard_scale` 回到 1.0（单元测试直接驱动 `Escape` 状态机）。
    status: pending

  - id: ac-006
    summary: 无意义的轴组合被具名拒绝，空死路梯诚实失败
    type: runtime
    pass_when: |
      `Rosenbluth { Enumerate }` 配 `Continuum` → `Grow::validate` 返回 `GrowError::AxisMismatch`
      并点名两条轴；`Escape` 列表为空的 CbmcGrow 在稠密算例上返回 `converged == false`
      且不 panic、不无限循环（有 max_loops 上限）。
    status: pending

  - id: ac-007
    summary: Python 镜像与文档同步
    type: code
    pass_when: |
      不增加 `WalkGrow` pyclass。`CbmcGrow` 无 `with_relax / with_serial /
      with_void_bias`，有 `with_schedule / with_reseed_void_bias / with_initial_hard_scale`；
      `LatticeGrow` 有 `with_schedule`；`.pyi` / `_protocols.py` 同步；
      docs/python/guide/growth.md 不含 `WalkGrow` 一节；
      CLAUDE.md 与 docs/architecture.md 的 `src/grow/` 行按新目录更新；`cargo doc` 零警告。
    status: pending
