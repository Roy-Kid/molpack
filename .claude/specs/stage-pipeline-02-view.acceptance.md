---
slug: stage-pipeline-02-view
criteria:
  - id: ac-001
    summary: RigidView owns the rigid DOF and nothing else; PlacementsMut is gone
    type: code
    pass_when: |
      `src/context/rigid_view.rs` 定义 `pub struct RigidView`（`Debug + Clone`），字段恰为
      `x` / `nmol`，方法恰含 `fresh` / `nmol` / `com` / `set_com` / `euler` / `set_euler` /
      `as_slice` / `as_mut_slice` / `write_xcart` / `install_seed` / `capture_from_xcart`，
      文件 ≤ 300 行；`grep -rn 'PlacementsMut' src/ python/src/ docs/` 无命中。
    status: pending

  - id: ac-002
    summary: context never names entry; the seed arrives as plain data
    type: code
    pass_when: |
      `grep -n 'use crate::entry' src/context/rigid_view.rs` 无命中，
      `grep -rn 'use crate::entry' src/context/` 无命中；
      `grep -n 'use crate::' src/context/rigid_view.rs` 不含 `gencan` / `grow` / `initial`；
      `install_seed` 的签名是 `(x: &[F], coor: &[[F; 3]], ctx: &mut PackContext) -> Self`，
      `Placements` 的拆包发生在 `src/gencan/entry.rs`。
    status: pending

  - id: ac-003
    summary: xcart is the single home before anything captures from it
    type: code
    pass_when: |
      `src/grow/driver.rs` 的 abort 补齐循环（今天的 `:352-357`）内含把
      `chain.coords` 写入 `sys.xcart` 的同步，与 `src/grow/lattice/mod.rs:262` 同形；
      `RigidView::capture_from_xcart` 只从 `ctx.xcart` 读坐标
      （`grep -n 'chain.coords\|d.coords' src/context/rigid_view.rs` 无命中）。
    status: pending

  - id: ac-004
    summary: The abort writeback is bitwise unchanged by the xcart move
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- grow_abort_writeback_golden` 通过：
      `EarlyStopHandler` 触发 abort 后的 `PackResult::positions()` 与测试内硬编码字面量逐位
      相同；该金标捕获自本 spec 之前的构建（注释记录），并在 xcart 同步与
      `capture_from_xcart` 落地之后仍绿；`grow_abort_keeps_bonded_geometry` 断言未改且全绿。
    status: pending

  - id: ac-005
    summary: One writeback derivation, one xcart rebuild
    type: code
    pass_when: |
      `grep -rn 'init_xcart_from_x' src/` 无命中；
      `grep -n 'use crate::initial' src/entry/mod.rs` 无命中；
      `grep -rn 'set_euler' src/grow/` 无命中（两处写回块已被 `capture_from_xcart` 取代）；
      `grep -n 'seeded' src/context/rigid_view.rs` 只命中说明该标记为何**不在**视图上的
      rustdoc 散文，不命中任何字段或方法定义。
    status: pending

  - id: ac-006
    summary: Seeded continuity and the three entries stay bitwise
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests` 通过（既有 RED
      `grow_cg_kremer_grest_c_inf` 除外，未被修改或跳过），`tests/packer.rs::seeded_run_contract`
      与 `free_chain_push_off_deterministic`、`tests/grow.rs` 的 `grow_deterministic_same_seed` /
      `grow_copy_stream_independence` / `lattice_grow_bead_chain_constructive` 断言未改且全绿；
      `cargo test -p molcrafts-molpack --release --features io --test examples_batch -- --ignored`
      五例通过。
    status: pending

  - id: ac-007
    summary: Regression scenario reproduces hard-coded rigid-view goldens
    type: runtime
    pass_when: |
      `cargo test -p molcrafts-molpack --lib --tests -- rigid_view_regression_xcart_and_capture_golden`
      通过：已知 `(com, euler, coor)` 经 `write_xcart` 得到的坐标与测试内硬编码字面量在
      1e-12 内相等，反向 `capture_from_xcart` 回到同一 `(com, coor)`；无第三方运行时。
    status: pending
---

# Acceptance criteria

- **ac-001 / ac-002** 是边界门：一个类型持有刚体自由度（且只有它），`context` 不认识 `entry`。
- **ac-003 / ac-004** 把"`xcart` 先成为唯一的家"变成可查的结构 + 数值证据。
- **ac-005** 是"一个事实一个家"的结构门：一个写回推导、一个 xcart 重建，视图上没有第二份
  "已放置"表示。
- **ac-006** 是数值门：种子接续与三个入口逐位不变。
- **ac-007** 是本 spec 的回归场景。
