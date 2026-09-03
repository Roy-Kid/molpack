---
slug: stage-pipeline-03-state
criteria:
  - id: ac-001
    summary: PackState wraps PackContext without moving a field, and stays crate-private
    type: code
    pass_when: |
      `src/context/pack_state.rs` 定义 `pub(crate) enum Placed { None, All }` 与
      `pub(crate) struct PackState`，字段恰为 `ctx` / `placed` / `rigid`（无 `topology`，
      `rigid` 非 `Option`），文件 ≤ 300 行；
      `grep -n 'PackState\|Placed' src/lib.rs` 无命中（本 spec 不加公开面）；
      `src/context/pack_context.rs` 的字段集合与行数相对本 spec 之前未增；
      `grep -n 'use crate::' src/context/pack_state.rs` 不含 `gencan` / `grow` / `entry` /
      `initial` / `refine`。
    status: pending

  - id: ac-002
    summary: One home and one idiom for the unscaled verdict
    type: code
    pass_when: |
      `grep -rn 'fn evaluate_unscaled' src/` 只命中 `src/context/pack_state.rs`；
      `grep -rn '\.scale2 = ' src/` 的命中恰为三条——`src/initial.rs:330` 与
      `src/context/pack_state.rs` 内该 helper 的置位与对称还原两行——`src/grow/` 与
      `src/gencan/` 下零命中（`src/context/pack_context.rs:370-371` 是结构体字面量默认值，
      不匹配该模式；`src/objective.rs:414,622` 是 `let scale2 = sys.scale2;` 读取，同样不匹配）；
      `src/gencan/phases.rs` 不再定义该函数，只调用它。
    status: pending

  - id: ac-003
    summary: The helper restores scale/scale2 as symmetrically as it restores radius
    type: runtime
    pass_when: |
      `tests/context_pack_state.rs` 中一条测试在人为把 `ctx.scale` / `ctx.scale2` 设为非默认值
      后调用 `evaluate_unscaled`，断言返回后两个字段逐位等于调用前的值，且 `ctx.radius`
      逐位还原；该测试在本 spec 之前的实现下不可能通过（RED 记录在提交信息里）。
    status: pending

  - id: ac-004
    summary: Value parity on both paths
    type: runtime
    pass_when: |
      `tests/context_pack_state.rs` 的两条值等价测试通过：GENCAN 夹具（`radius` 已被 discale
      缩放）与生长夹具（`radius == radius_ini`）在合并前后返回的 `(f_total, fdist, frest)`
      与调用后的 `sys.radius` 内容逐位相同（`to_bits()` 断言）。
    status: pending

  - id: ac-005
    summary: No second placed bitset
    type: code
    pass_when: |
      `grep -n 'Vec<bool>' src/context/pack_state.rs` 无命中；`PackState` 不含 per-atom
      放置数组，只含 `Placed` 形状标记。
    status: pending

  - id: ac-006
    summary: The module doc states why the swap is safe and books the surface removal
    type: docs
    pass_when: |
      `src/context/pack_state.rs` 的模块 / 函数文档写明：`scale` / `scale2` 的三个写点
      （含 `src/initial.rs:329-330`）都写默认值；两组字段对称存还；`fdist` 由 `radius_ini`
      无条件计算（`src/objective.rs:292-300`）而 `frest` 读 `scale` / `scale2` 但不读半径，
      故半径互换只移动 `f_total`；并明账 `molpack::gencan::phases::evaluate_unscaled` 这条
      已发布路径被撤下、文档站两处（`docs/architecture.md:37`、`docs/extending.md:460`）由
      07 跟随；`cargo doc -p molcrafts-molpack --no-deps` 零警告。
    status: pending

  - id: ac-007
    summary: The suite and the Packmol regression stay green
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` 通过（既有 RED
      `grow_cg_kremer_grest_c_inf` 除外，未被修改或跳过）；
      `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored`
      五例通过。
    status: pending

  - id: ac-008
    summary: Regression scenario reproduces hard-coded verdict goldens
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- pack_state_regression_unscaled_verdict_golden`
      通过：固定夹具（6 二聚体 / 20 Å 盒 / seed 7）的 `(f_total, fdist, frest)` 与测试内硬编码
      字面量在 1e-12 内相等；注释记录金标捕获自本 spec 之前的构建；无第三方运行时。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-005** 钉住"包裹不抽取"、"无 `topology` 空字段"、"不造第二份放置位集"、"不加公开面"。
- **ac-002 / ac-003 / ac-004** 是未缩放裁决合并的结构门、对称性门与数值门。
- **ac-006** 让金标断言建立在被写下来的等价性论证上，并把公开路径的撤下记成明账。
- **ac-007 / ac-008** 是回归与本 spec 的回归场景。
