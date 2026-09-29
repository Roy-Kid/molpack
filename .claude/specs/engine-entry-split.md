# engine-entry-split — 按算法拆分引擎入口

状态:DRAFT(2026-09-01)。Owner 指令:删除大一统 `Molpack`;pack 与 grow、以及不同 grow
各自独立入口;所有入口共享一个管理生命周期/handler 的 trait;体系可能带支链与不同拓扑;
不保留历史包袱。architect design-mode 审查已完成(1C/9H/8M,两个硬门槛,见 Design)。

> **2026-09-29 名称对照**:本 spec 已 DONE(2026-09-01),正文保留当时的设计名。现行名:结果类型 → 冻结的
> `State`(`src/entry/result.rs`);跨入口接续 → `GenCanPack::with_restart(&State)`;`Solver` / `GencanSolver`
> → `Stage` / `GencanStage`(`src/stage.rs`、`src/gencan/solver.rs`);入口 trait 的 `solver()` 钩子 → `StageFactory::stages`;
> `examples_batch` harness 于 2026-09-20 删除,GENCAN 回归改由五个 `--example pack_<name>` 程序与
> `gencan/entry.rs::gencan_entry_is_deterministic`、`gencan/solver.rs::gencan_solves_a_small_pack_on_the_seam` 承担。

**Amends `NOTES.md 2026-08-28 §3`**:方法选择从 target 级(`Target::with_method`)迁至入口级;
"绝不静默降级到另一方法"条款原样存续(入口即方法,错配模板 named-reject)。顺带解决
`driver.rs:171` 的 per-target 旗标被 `any()` 折叠问题(选择上移后不复存在)。

## Goal

用户面从一个 26 旋钮的 `Molpack` 变成按算法命名的入口类型:**`GenCanPack`**(刚体 GENCAN)、
**`CbmcGrow`**(现连续构象偏置生长)、**`LatticeGrow`**(金刚石格相生长,见
lattice-growth-phase spec)。所有入口实现同一个生命周期 trait:trait 拥有
①校验 ②盒子解析(density/pbc/cell)⑤组装/报告 与全部横切关注(handlers 括号、seed、
tolerance、budget、log);入口只提供一个抽象钩子——返回它的 `Solver`。树状(支链)模板是
两个 grow 入口的一等输入;环模板显式拒绝。`Molpack`、`PackMethod`、`compose.rs`、Python
`_engines.py` 全部删除,无别名、无过渡层。

## Design(architect 裁定,两个硬门槛)

**门槛 1 — `GencanSolver: Solver` 必须先落地。** `Solver` 缝今天只有一个实现者
(`GrowthSolver`);GENCAN 实际是 `Molpack::run_gencan_stages` 的 12 参内在方法——"对等"
只存在于文档。第一个 RED:`src/gencan/solver.rs` 的 `GencanSolver`(12 参收编为字段),
`examples_batch` 保持绿。此后入口 trait **不是第二道缝**,而是生命周期模板:

```rust
pub trait PackEngine {  // 名称待定夺;终态动词唯一
    fn settings(&self) -> &PackSettings;
    fn solver(&self, targets: &[Target]) -> Result<Box<dyn Solver>, PackError>;
    fn run(self, targets: &[Target]) -> Result<State, PackError> { /* provided */ }
}
```

- `run(self, …)` **按值消费**——`pack(&mut self)` 里 `mem::take(self.handlers)` 导致二次调用
  静默无头运行的潜伏 bug,由类型系统根除(one-shot 类型强制)。
- 终态方法唯一,返回 `State`(`frame` 是其字段);不复刻 `pack`/`pack_with_report` 双门面。
- handler 存储、on_start/on_finish 括号、should_stop 抽取、log-handler 注入全部住在 trait 的
  provided `run()` 里(单一 `LogSpec` 驱动);`&mut [Box<dyn Handler>]` 照旧流入
  `Solver::solve` 做逐步回调。`on_phase_start`/`on_phase_end` 文档标注为
  GENCAN-only 可选钩(`on_inner_iter` 已于 2026-09-29 删除);grow 入口不得引入平行的 stage 观察者。

**门槛 2 — push-off 链式必须显式化。** 今天 grow 不收敛时 packer 内部悄悄接
`run_gencan_stages(push_off=true)`——单方法入口若保留即隐藏第二算法,若丢弃即熔体密度
静默倒退。裁定:`CbmcGrow`/`LatticeGrow` 如实报告 `softened`/`converged`;**用户显式链式**:

- push-off 链式:同一批 target 保持 **free**,喂 `GenCanPack`(推开分子间接触);
- 固定基质链式(原 compose 语义):`Target::fixed_from(&State)` 新原语(`src/target.rs`),
  长成的链变 fixed target 喂下一入口。两种链式各钉一个测试。
- `compose.rs` 删除时三项自有行为的去处:合并判据无需新家(stage-B 上下文含 fixed 链,
  其 fdist/frest 已裁决全系统;`compose.rs:180` 的 `max()` 是赘余,记录之);盒子经由共享
  `PackSettings` 值流入两个入口(不许变成"用户传两遍");输出顺序 = 用户列表序。
  **不引入 `Pipeline` 类型**(无调用方,earn-complexity)。

**`PackSettings`(one-home)**:tolerance、precision、discale、seed、box/density/cell、全局
restraint 广播——被共享基建消费的旋钮住一个值类型,trait 提供一次 `with_*` 转发。
GENCAN 专属旋钮(`inner_iterations`、`init_passes`、`init_box_half_size`、`perturb*`、
`avoid_overlap`)只在 `GenCanPack` 上,不得上 trait。

**拓扑**:树状模板两个 grow 入口一等支持(`InternalTree` 已树分解,
`internal_roundtrip_branched` 钉住);**环模板今天被静默摊平**(`bfs_order` 丢弃闭环键,
开链生长,无错误)——本 spec 为**两个** grow 入口加 `GrowError::RingTemplate` 拒绝。

## Non-goals

- 不做 `Pipeline`/工作流类型;链式是两行用户代码 + 一个 `fixed_from` 原语。
- `.inp` 不新增任何生长关键字——**`.inp` 按构造即 GENCAN-only**,记录为不变量;未来若开放
  必须是 script 级入口选择关键字,绝不 per-structure。
- 不保留 `Molpack`/`Packer`/`Grower` 任何别名或废弃壳(owner:不留历史包袱)。
- 不动 `ff` 边界(四原则 §1)。

## Public surface

**Rust 新增**:`PackEngine` trait + `PackSettings`(`src/entry.rs`);`GenCanPack`
(`src/gencan/entry.rs`)、`CbmcGrow`(`src/grow/entry.rs`)、`LatticeGrow`
(`src/grow/lattice/entry.rs`);`Target::fixed_from(&State)`;`GrowError::RingTemplate`;
`MolpackLogLevel` → `LogLevel`。
**Rust 删除**:`Molpack`、`pack`/`pack_with_report`、`PackMethod`、`compose.rs`、`packer.rs`。
**Python**:删除 `_engines.py`(`Packer.on()`/`Grower.run()` 两动词一操作的手写糖);入口
pyclass 与 Rust 同名 1:1,单一终态动词;共享旋钮绑定一次(trait 级 py 基座),不三处重复;
`_protocols.py`/`.pyi` 重生成。破坏面清单:`python/src/packer.rs`(507 行重写)、
`GrowConfig` 旋钮迁 `CbmcGrow`、`PackMethod` py-enum 删除、7 个测试文件 + 6 个示例迁移;
`with_progress` 仅存于 Python 侧,迁 `LogSpec`。
**CLI**:`ScriptPlan.packer`/`BuildResult.packer` 字段改 `entry: GenCanPack`
(`build.rs` 只设 tolerance/avoid_overlap/seed/pbc/cell,证实 .inp 只降到 GENCAN)。

## Module placement(packer.rs 1901 行 = 2.4× 硬上限,拆分是本 spec 的一部分)

| 新位置 | 内容 | 预算 |
|---|---|---|
| `src/entry.rs` | `PackEngine` + `PackSettings` + `State` | ≤300 |
| `src/entry/setup.rs` | density/pbc/cell 解析(原 packer.rs:533-650)+ restraint 广播 | ≤300 |
| `src/context/build.rs` | `build_context`(原 836-1047) | ≤250 |
| `src/gencan/solver.rs` | `GencanSolver: Solver`(原 run_gencan_stages) | — |
| `src/gencan/phases.rs` | run_phase / run_iteration / evaluate_unscaled(原 1512-1869) | — |
| `src/gencan/entry.rs`、`src/grow/entry.rs`、`src/grow/lattice/entry.rs` | 三个入口 | 各≤300 |

删除:`src/packer.rs`、`src/compose.rs`。

## Numerical contract

- GENCAN 数值零变化:当时的 `examples_batch`(release, --ignored)全绿是门槛 1 的验收(该 harness 已删,现行守门见文首对照);
  重构是代码搬运,目标函数/优化器路径逐位不动。
- grow 确定性测试(同种子逐位、copy 流独立、mixed 逐位)在新入口下原样通过
  (混合逐位测试改写为显式 `fixed_from` 链式的等价断言)。
- push-off 链式的行为等价:`CbmcGrow` 不收敛样例 + 显式 `GenCanPack` free 链式,产物与
  今日隐式路径逐位一致(迁移期一次性对照,之后旧路径删除)。

## Test plan

- RED-1:`GencanSolver: Solver` 单测(现 `gencan/solver.rs::gencan_solves_a_small_pack_on_the_seam`)+ 当时的 `examples_batch` 绿(先于一切入口代码)。
- 入口生命周期:consume-by-value(二次 run 编译不过,doc-test 展示);handler 括号/
  should_stop/log 注入在三个入口上各一冒烟(共享 provided run,测一处逻辑三处接线)。
- 链式:push-off 链式与 fixed-matrix 链式各一集成测试(后者吸收
  `test_fixed_only.py::test_fixed_matrix_plus_grow_target` 语义)。
- `RingTemplate` 拒绝:两个 grow 入口 × 环模板(今天静默摊平的用例转为 RED)。
- 分支模板:`branched_parts()` 几何过 `CbmcGrow` 全链路(现仅内坐标往返有测试)。
- Python:入口 1:1 冒烟 + 迁移后的 7 测试文件全绿。

## Doc plan

- `docs/architecture.md`:模块表按新布局重写,**并纠正当时文档里的陈旧项**(文档仍列出早已不存在的 restraint.rs、relaxer.rs、
  cell.rs、api/ 模块;83-86 行的虚假对等声明删除);CLAUDE.md 模块表同步;
  `docs/python/` 全部入口示例迁移。

## Risks / open questions

1. 一次性破坏面大(Rust+Python+CLI+7 测试+6 示例)——顺序上 RED-1 → trait+GenCanPack →
   CbmcGrow → 删除旧面 → LatticeGrow(lattice spec 依赖本 spec 落地)。
2. trait 名与终态动词名(`PackEngine::run`?)——owner 定夺;"GenCan" 是优化器本名
  (Birgin & Martínez),不触 no-"packmol" 规则。
3. `short_tolerance` 归 `PackSettings` 还是 GENCAN 专属——按消费方(build_context)裁定,落地时核对。
4. lattice-growth-phase.md 的 `PackMethod::GrowLattice` 过渡挂载点作废,该 spec 已同步改为
   `LatticeGrow` 入口。

## 交叉引用

- 依赖方:`lattice-growth-phase.md`(LatticeGrow 入口在此 trait 上)。
- 修订:`NOTES.md 2026-08-28 §3`(见文首 Amends)。

## 落地修正案(2026-09-01,实现记录;门槛 2 已按原裁定收尾)

1. **门槛 2:free-target 链式已落地**(`placement-seeding.md`):
   `State` 携带放置解(`Placements`:x 向量 + 逐副本居中构象 + cell),
   `GenCanPack::with_restart(&grown)` 逐位原样接续——同批 free target、
   `initial()` 跳过、movebad 关闭、盒子随 seed 流动(不传两遍)。与过渡形态
   `with_push_off` 的逐位等价由一次性迁移测试
   `free_chain_matches_push_off_bitwise` 证明后,`with_push_off` 与
   `PackEngine::after_solve` 钩子一并删除(过渡形态自身此前已逐位等价于
   更早的隐式路径——等价链完整)。行为不变式常驻:
   `free_chain_push_off_starts_from_grown_state`(交接逐位)、
   `free_chain_push_off_deterministic`、`seeded_run_contract`。
2. **`after_solve` 钩子:已删除**(唯一实现者随门槛 2 收尾退役)。
3. **GENCAN 缺省一处安家**:`GencanSettings::default()`;`push_off` 并入
   `GencanSettings`(非 `GenCanPack` builder 旋钮,由链式代码设置)。
   `GencanSolver::new` 不再接收 optimizer 参数——optimizer 只经
   `GenCanPack::with_optimizer` 接入 `src/gencan/`,验收门
   `grep -rn 'molrs::optimize\|optimizer::' src/grow/` 保持零命中
   (2026-09-29 起 optimizer 接缝不再受 `ff` 门控,`ff` 只透传 `molrs/ff`,
   原 `cfg(feature = "ff")` 判据失效)。
4. **Python 破坏面补记**:`with_lammps_output`、`with_xyz_output` 随
   `Molpack` 删除(前者 ≡ `with_log_level("progress")`,后者由
   `with_handler` 覆盖);`with_log_level(str)`/`with_log_frequency`
   保留在两个入口上。共享 builder 由 `entry_pymethods!` 宏一次绑定,
   两入口不再手抄。`with_cell`/`with_short_tolerance` 尚未上 Python
   (旧绑定也没有;需要时随用例补)。
5. **`LogSpec`** 退出 crate 根重导出(无公开签名使用),经
   `entry::LogSpec` 可及。
