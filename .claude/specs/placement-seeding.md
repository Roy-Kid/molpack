# placement-seeding — 跨入口延续放置解,门槛 2 按原裁定落地

状态:LANDED 2026-09-01(一次性逐位对照 `free_chain_matches_push_off_bitwise`
通过后 `with_push_off`/`after_solve` 已删除;常驻测试见 Test plan)。
现行公开名(2026-09-04 `result-as-state`):结果叫 `State`,接续入口叫
`GenCanPack::with_restart`;下文已按现行名改写。

## Goal

一个入口的运行结果能被下一个入口**逐位原样**接续:`CbmcGrow` 如实报告
不收敛后,用户把**同一批 target(保持 free)**喂给
`GenCanPack::with_restart(&grown)`,GENCAN 相位在长成的坐标上做刚体推开
(Auhl slow push-off)。`CbmcGrow::with_push_off` 与
`PackEngine::after_solve` 钩子随之删除——增长入口内不再藏第二算法。

## 核心观察(为什么不用 frame 重建)

从装配后的 frame 反推 (coor, x) 需要重算 COM:
`(p − com) + com ≠ p`,逐位等价即失。放置解必须**原样携带**——
`State` 在 Stage ⑤ `init_xcart_from_x` 之后捕获 `(x, coor_free, cell)`
的逐位快照,seeded 运行原样注入(zero-conversion chaining,与
chain-growth-solver Design §2 同一原则)。

## Public surface

Rust:

- `State` 新增 `pub(crate) placements: Placements`(crate 私有字段;
  对外不可见,不改公开构造面——`State` 本就只由 `run()` 构造)。
- `entry/result.rs`:`pub(crate) struct Placements { x: Vec<F>,
  coor: Vec<[F;3]>, copy_atoms: Vec<usize>, cell: SimBox }` —— free 副本的
  (COM|Euler) 打包向量、逐副本居中参考构象(xcart 序)、逐副本原子数
  (校验指纹)、本次运行安装的 simbox。
- `GenCanPack::with_restart(&State) -> Self`(GENCAN 专属,不上 trait
  ——增长入口自己造初态,seed 对它无意义)。同时把 seed 的 cell 写入
  `settings.cell = CellDecl::Matrix{…}`:盒子随 seed 流动,**用户不传两遍**
  (门槛 2 判词);用户若又声明 box/density/cell,由**既有**互斥错误具名
  拒绝(不新增静默优先级)。
- `PackEngine::prepare` 签名扩为 `(…, x: &mut [F], …)` —— 入口级上下文
  准备天然包含放置种子;`CbmcGrow::prepare` 忽略 x。
- 新错误:`PackError::SeedMismatch { expected, got }`(seed 的 free 原子
  形状与本次 targets 不符)。
- 删除:`CbmcGrow::with_push_off`、`PackEngine::after_solve`、
  (随之)`grow/entry.rs` 对 `gencan` 的全部引用。

Python 镜像:`GenCanPack.with_restart(state)`;`CbmcGrow.with_push_off`
删除;`.pyi` 同步。

## Semantics

1. **注入点**:`GenCanPack::prepare`(seed 在手时)——
   `install_simbox_and_grid(sys, seed.cell, radmax, discale, ntotat_free)`
   (与 `CbmcGrow::prepare` 同一套安装),然后
   `sys.coor[..ntotat_free] ← seed.coor`、`x ← seed.x` 原样拷贝。
2. **求解模式**:seed 在手 ⇒ `GencanSettings.push_off = true`:跳过
   `initial()`(它会重掷全部 COM/Euler,把长成的链传送走),经
   `init_xcart_from_x` 物化 xcart(与上一入口 Stage ⑤ 同一算式、同一
   输入 ⇒ 逐位同一 xcart),movebad 关闭——分子只被刚体下降推开。
3. **额外 fixed target 允许**:seed 只覆盖 free 块;free 形状必须与
   seed 逐副本一致(`SeedMismatch` 否则),其后可以追加 fixed 基质。
4. **verdict 诚实**:链上每段自报 `converged/fdist/frest/softened`;
   增长段的 `softened` 留在它自己的 `State` 上,推开段 GENCAN 恒 0
   ——比旧隐式路径把两段搅在一起更诚实。
5. **frame 盒子**:`run()` 的 simbox 盖章规则扩为:`space.pbc` 优先,
   否则**声明过 cell**(含 seeded 注入)时盖 `space.cell`;仅约束推断的
   盒子(用户没声明)照旧不盖。

## Numerical contract

- 迁移期一次性逐位对照(spec 原句):`CbmcGrow` 不收敛样例上,
  `with_restart` free 链式(同 tolerance/discale、GenCanPack 的 seed =
  增长入口的 seed、同 max_loops)产物与 `with_push_off(true)` **逐位一致**
  ——证毕后删除 `with_push_off`,对照测试转为行为不变式(交接连续性
  bitwise、刚体性、verdict 诚实、同种子确定性)。
- 逐位成立的机理链:seed 原样拷贝 ⇒ 同 (coor,x);SimBox 由同一 (h,origin,
  pbc) 重建 ⇒ 同网格;同 targets+tolerance ⇒ 同 radii/maxmove;GencanSolver
  自持 RNG(seed 相同)⇒ 同随机流。

## Test plan

- RED:`free_chain_matches_push_off_bitwise`(当时的集成测试,删除前的
  一次性对照)。
- 迁移:`grow_push_off_starts_from_grown_state` → 链式版(HandoffProbe 在
  第二段 `on_initialized` 捕的 xcart 与 `grown.positions()` 逐位相等);
  `grow_push_off_deterministic` → 链式版。(二者与 `tests/` 一并于 2026-09-20
  删除;逐位交接比较不再常驻。)
- `SeedMismatch` 具名拒绝;seed + `with_periodic_box` 撞既有互斥错误;
  seeded + 追加 fixed 基质一例(§3)——现由
  `src/grow/tests/entry.rs::seeded_run_contract` 一并钉住。
- Python:`with_restart` 链式冒烟 + one-shot 语义不变;7 文件全绿。

## Out of scope

- 增长入口接收 seed(无意义);LatticeGrow 的 seed 语义随其 spec。
- 跨进程持久化 placements(mrec 序列化)——需要时另立。
