---
slug: chain-growth-solver
criteria:
  - id: ac-001
    summary: src/grow/ 与 src/stage.rs 不依赖 ff feature，先验是纯几何数据
    type: code
    pass_when: |
      `cargo build -p molcrafts-molpack`（default feature，不带 io/ff/cli/rayon）
      成功，且 src/grow/（含 prior.rs）与 src/stage.rs 均已编入。
      `grep -rn 'molrs::ff\|molrs::optimize' src/grow/ src/stage.rs`
      无命中（接缝文件现为 `src/stage.rs`；crate 内已无
      `cfg(feature = "ff")`，`ff` 只透传 `molrs/ff`，故不再以该 cfg 为判据）。TorsionPrior / AnglePrior 的构造只收数据（角度、权重、C∞、
      持久长度），不收任何力场对象。
      `cd python && maturin develop --release` 后 `import molpack` 能取到生长入口，
      wheel 的 feature 集合未新增。
    status: pending

  - id: ac-002
    summary: 生长是 Solver 接缝上的平级算法，方法选择是 per-target 的
    type: code
    pass_when: |
      （2026-09-29 按现行接缝改写：`Solver`/`PlacementsMut`/`SolveOutcome`、
      `Target::with_method`/`PackMethod`、`pack_with_report` 已随 engine-entry-split
      与 stage-pipeline 删除。）
      src/stage.rs 定义 `pub trait Stage`（`name`/`requires`/`guarantees`/`run`，
      `run` 收 `&mut PackState`，返回 `Result<StageOutcome, PackError>`；
      Budget / StageOutcome 标 `#[non_exhaustive]`）。方法选择由调用方挑入口
      （`GencanPack` / `CbmcGrow` / `LatticeGrow`，或 `Pipeline::with_stage` 组合），
      不存在全局 solver 开关；`GrowConfig` /
      `GrowError` 定义在叶子文件 src/grow/config.rs（`grep -nE 'use crate::(target|entry|context)' src/grow/config.rs` 无命中）；
      公开结果是冻结的 `State`（`frame`/`fdist`/`intra`/`frest`/`converged`/`degraded`），
      原 `softened` 由 `State::degraded` 承载（各阶段 `StageOutcome::softened` 之和）。
      `grep -rn 'pgencan\|run_phase\|run_iteration' src/grow/ src/stage.rs`
      无命中——solver 永不调用 pack 内部代码，生长与 GENCAN 并列不嵌套。
      `Method::Grow` 配给原子数 < 3 或无键模板返回具名错误，不静默退化为刚体。
    status: pending

  - id: ac-003
    summary: 内坐标往返精确，且任意自由变量下成键几何不变（环键误判探测器）
    type: runtime
    pass_when: |
      src/grow/tests/internal.rs 中（`internal_roundtrip_{linear,branched,ring}`、
      `internal_random_vars_preserve_bonded_geometry`、`internal_rigid_molecule_no_vars`）：
      (i) 对模板帧构建 InternalTree 后用模板扭转值重建，
      与模板坐标逐点比较 ‖Δ‖∞ < 1e-9；(ii) 用**随机**自由变量重建后，所有
      键长、键角、非自由二面角、以及环闭合键长仍等于模板（1e-9 / 1e-9 rad）。
      (ii) 是 (i) 抓不住的两类 bug（环内键被误判为自由变量、步分组错误）的
      唯一机械探测器。三个边界用例都过：带支链、含环、无可旋转键分子。
    status: pending

  - id: ac-004
    summary: 无重叠是构造保证，且报告数字来自共享 objective
    type: runtime
    pass_when: |
      生长产物的 `State::frest == 0.0`（restraint 无违反；原
      `validation::validate_from_targets` 已于 2026-09-29 删除，裁决读 `State`）；
      `State::fdist == 0.0`（严格等零）；`State::degraded == 0`
      （公开载体即此字段，原 `softened`）；
      独立复算的最小非排除原子间距（最小镜像下）≥ tolerance。
      SolveOutcome 的 fdist/frest 由收尾时对共享 objective（Constraints 入口）
      的一次评估产出，src/grow/ 中不存在自行赋值最终 fdist/frest 的路径。
      至少在 dp=100×25 / ρ=1.0 与 dp=200×25 / ρ=1.0 两个算例上成立——
      这两个正是刚体路径实测停在 fdist=3.3252（最近对 0.82 Å）与
      60 loops 不收敛的地方。
    status: pending

  - id: ac-005
    summary: 密度归 Molpack 层，直接命中，无默认值
    type: runtime
    pass_when: |
      `Molpack::with_density(rho)` 存在（GrowConfig 上没有 density 参数）；
      盒长由**全部** targets（含混合体系的 Gencan targets）的总质量解析算出；
      质量默认查元素、`Target::with_mass` 可覆盖；元素为 "X" 且无覆盖 →
      具名错误不猜；
      装完后由输出帧 box 体积与总质量算出的实际质量密度与 rho 相对偏差 < 1e-6。
      与 with_periodic_box / with_cell 同时给 → 报错；含 Grow target 而
      既无盒子也无密度 → 报错，不回退到任何隐含值。
    status: pending

  - id: ac-006
    summary: 链统计阶梯——先验从弱到强逐级可测，熔体 Rg 达标且每链各异
    type: scientific
    pass_when: |
      (1) 孤立链、TorsionPrior::Uniform、无排除：C∞ = 2.00 ± 0.1
      （自由旋转链解析值，NeRF + 采样器的回归基线）。
      (2) 孤立链、RIS 先验 p_t = 0.645（由 three_state_from_c_inf(5.5, 四面体角)
      解出）：C∞ = 5.5 ± 0.3。
      (3) 熔体 dp=200（M=8829）× 25、rho=1.0、RIS 校准先验：平均 Rg 相对
      Flory 值 34.4 Å（<R²>₀/M = 0.805 Å²·mol/g）偏差 ≤ 10%，且严格优于
      模板构象的 +34%；各 copy Rg 互不相同（max − min > 1 Å）——对照实测病症
      min = mean = max = 23.1 Å。
      (4) 链序无偏：生长编号前 5 条与后 5 条的平均 Rg 差 < 5%。
      Rg 一律键图展开（bond-unwrapped）后计算。
    status: pending

  - id: ac-007
    summary: 约束硬拒绝——frest 与 fdist 一样是构造保证
    type: runtime
    pass_when: |
      一个带 `InsideSphereRestraint` 的生长算例：全部原子落在球内，
      `State::frest == 0.0`（严格零，非 < precision——候选在 r.f > 0 时
      被拒绝，不是被惩罚）。
      `grep -rn 'trait .*Restraint' src/grow/` 无命中——生长不定义新约束 trait，
      只调用现有 `AtomRestraint::f`。
    status: pending

  - id: ac-008
    summary: grow 与 gencan 可零转换串联
    type: runtime
    pass_when: |
      生长产物（`sys.coor` + `x`）直接交给 GENCAN 路径继续，中间不存在任何坐标
      变换代码；串联后第一轮的 fdist 不高于串联前。
      写回契约可断言：每个 copy 的 `x` 欧拉分量为 (0,0,0)，COM 等于该 copy
      `coor` 的质心（1e-9）。
    status: pending

  - id: ac-009
    summary: 既有 Packmol 路径逐条不变
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --features cli,ff,rayon` 全绿；
      五个官方 Packmol 例子（`cargo run --release --features io --example pack_<name>`）
      在固定 seed 下收敛且 `State::frest == 0`（原 `examples_batch` harness 与
      validation 已删除）。
      在 Task 1（只立接缝 + per-target dispatch）之后与全部任务完成之后各跑一次，
      两次结果一致。接缝是新增分支，不是对既有算法的改写——这是守门条件。
    status: pending

  - id: ac-010
    summary: dp=200×25 / rho=1.0 单线程跑完，四个数字入库
    type: performance
    pass_when: |
      35050 原子（dp=200，1402 原子/链，25 条），rho=1.0（L=71.6 Å），
      单线程完成且不返回 DeadEnd，耗时 ≤ 120 s。
      实测耗时、retractions、relax_regrowths、softened 四个数字写进对应
      commit message。对照基线记在同一处：刚体路径该体系 rho=0.5 下 60 loops
      不收敛（fdist 停在 2.97）。
      门限 120 s 是工程上限：复杂度 O(N_atoms × n_trials)，真实数字预期低一个
      量级；实测逼近上限说明回撤率异常，属设计问题而非调参问题。
      （注：ρ=1.0 是工程测试点；PEO 实验熔体密度 1.060 g/cm³ @ 353 K。
      用户裁定 2026-08-28：不追求精准密度——产物后续走 MD NPT 弛豫，
      与实验值 ~20% 内的偏差可接受；本条的 1e-6 判据属 ac-005，
      核验的是"命中用户声明的 ρ"这个算术，与实验值无关。）
    status: pending

  - id: ac-011
    summary: CG 合成链在 default feature 下达到目标持久统计
    type: scientific
    pass_when: |
      STRUCK 2026-09-29：KG 珠链 c∞ = 1.76 夹具随集成测试层于 2026-09-20 删除，无替代；
      采样链的 C∞ 契约现由 src/grow/tests/prior.rs 的
      `prior_uniform_freely_rotating_c_inf` / `prior_ris_calibrated_c_inf` 钉住。以下为历史文本。
      （原）程序化合成的 Kremer–Grest 型珠链模板（无 io、无 ff）：
      `AnglePrior::Wlc`（或 States）校准后，孤立链 c∞ = 1.76 ± 10%
      （Auhl et al. 2003 的 NRRW 目标值）。
      `exclusion_depth` 为 per-target 参数且该用例显式设置（不依赖 AA 默认 3）。
      角度成为采样自由度的路径（AnglePrior != Template）被该用例实际走到。
    status: pending

  - id: ac-012
    summary: 确定性与 RNG 契约——并行预埋在 v1 串行版即生效
    type: runtime
    pass_when: |
      同 seed 两次运行输出逐位一致。
      RNG 流按 (seed, copy, step) 哈希为独立流：单元测试断言从体系中移除一条链
      不改变另一条链的候选提议序列（轮快照语义 + 独立流的联合探测器）。
      `OverlapField::probe` 签名保持 `&self`。
      轮快照语义（轮内提议针对轮初快照、提交串行 + 增量复查）在
      src/grow/mod.rs 模块文档中写明为算法定义的一部分，不是实现细节。
    status: pending

  - id: ac-013
    summary: 混合体系——per-target 方法在同一次 pack 里共存
    type: runtime
    pass_when: |
      主线：PEO(`Method::Grow`) + 小分子(`Method::Gencan`) 算例——盒子按总质量
      定容，生长先行，链原子在刚体阶段作为 fixed 结构逐位不动，整体
      `State::frest == 0`（原 validation 已删）。
      Task 9 检查点缩容时的替代判据：混合方法返回具名错误（不是 panic、不是
      静默单方法），且后续 spec 已在 INDEX 立项。两种结果都算过，静默错配不算。
    status: pending
---

# chain-growth-solver — 验收

本文件是 `/mol:impl` 的完成契约。每条 `pass_when` 必须可被机械核验；
`type: scientific` 与 `type: performance` 由 `examples/pack_peo`（AA）产出的量化报告佐证
（原集成层的 CG 用例已于 2026-09-20 删除；先验统计见 src/grow/tests/prior.rs）。

判据对照基线来自两处：实测数据（Domain basis §一–四）与文献值
（Domain basis §五：C∞ = 2.00 解析值、PEO C∞ = 5.51、KG c∞ = 1.76、
Auhl 0.8σ push-off 边界）——不是估计值。

密度均匀性 E(d) 与内距曲线 ⟨R²(s)⟩/s 是**报告项**不是断言项（文献无可断言的
通用门限），随 pack_peo 输出，供人工审阅。

ac-006(3) 的偏差方向可预期：Theodorou–Suter 原文实测长程偏置生长使链相对
无扰尺寸收缩约 8–15%（expansion factor 0.85/0.92，spec §5.7），故熔体 Rg
预期从**下方**接近 34.4 Å；若实测超出 +10%，应怀疑先验校准或排除表，而不是
放宽容差。
