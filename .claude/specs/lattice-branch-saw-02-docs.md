---
title: lattice-branch-saw-02-docs — LatticeGrow 星形 PEO 的 Python 示例与文档
status: code-complete
created: 2026-09-05
grilled: true
chain: lattice-branch-saw (02 of 02; predecessor lattice-branch-saw-01-walk)
---

# lattice-branch-saw-02-docs — LatticeGrow 星形 PEO 的 Python 示例与文档

## Summary

在 `01-walk` 已使 `LatticeGrow` 接受树状（含分支）模板的前提下，把 Python 示例 `pack_peo_topo.pack_star` 改成调用方的一次显式选择：`make_star` → `Target`（`PEO_H_RADIUS`）→ `LatticeGrow` @ 2.0 Å → `GencanPack.with_restart` @ 2.0 Å；删除星形路径上的格相 try/except 探测和 `CbmcGrow` 主体，也不再打印 Auhl。文档只公布验收集（树含分支可生长、环为具名 `RingTemplate`、非四面体仍具名拒绝），不把拓扑写成算法分派，也不泄露 01 的 C3/BFS。`CbmcGrow` 仍是对等的树生长器；Auhl 减排斥体积配方留在既有 `CbmcGrow` 路径。本段在 01 code-complete 之前不可实施。

## Domain basis

本段不引入新方程。格上分叉 SAW 由前驱 `lattice-branch-saw-01-walk` 拥有。`LatticeGrow` 在占据守卫下生成；`with_tolerance` 是装饰后共享 objective 的尺子，不是连续 CBMC 的第一阶段减核。既有教师（`docs/python/guide/growth.md` 格相样例、`examples/pack_peo` 的 `lattice` 动词、`python/tests/test_grow.py::TestLatticeGrow.test_lattice_then_seeded_push_off`）一律用打包半径 2.0 Å 生成，再 `GencanPack.with_restart` @ 2.0 Å。把 `PEO_GROW_TOL`（默认 0.6 Å）抄到格相入口，会造出第二套更松的 `State` 裁决。

Auhl 减排斥体积生长 → 慢 push-off 只属于 `CbmcGrow` 路径（`test_auhl_two_stars_push_off_converges`、Rust `examples/pack_peo` `auhl`）。Auhl, Everaers, Grest, Kremer & Plimpton, *J. Chem. Phys.* **119**, 12718 (2003), doi:10.1063/1.1628670。

- 容差单位：Å。默认打包半径 / 装饰后尺子：2.0 Å。
- 密度单位：g/cm³。氢打包半径：`PEO_H_RADIUS` 默认 0.85 Å。

## Design

本段不改 Rust、不改 PyO3、不新增公开类型。公开入口仍是对等的 `LatticeGrow` / `CbmcGrow` / `GencanPack`。示例是一次 pick（P8），不是拓扑→算法分派，也不是求解器内部回退。树行走的所有权在 01；本段只翻转调用方与文档，使它们与 01 落地后的公开合同一致，而不是在 01 之前另写一套「已接受分支」的真理。

### pack_star

当前 `pack_star` 先用 `LatticeGrow` @ 2.0 Å 做 try/except 探测，失败后再走 `CbmcGrow` @ `PEO_GROW_TOL`（0.6 Å）并打印 `Auhl`。重写为唯一配方：

1. `polymer = make_star(dp, seed=seed)`（molrs SMILES + molpy `PolymerBuilder`，不手摆坐标）。
2. `target = _target(polymer, n_mol, "star-PEO")`（保留 `PEO_H_RADIUS`，默认 0.85 Å）。
3. `LatticeGrow(prior).with_seed(seed).with_tolerance(2.0).with_density(density)`（占据守卫默认 on；`max_loops=max(40, n_mol * 8)`）。
4. `GencanPack().with_restart(grown).with_seed(seed).with_tolerance(2.0)`（`max_loops=80`）。

删除：格相 try/except；`pack_star` 体内全部 `CbmcGrow`；任何 `print("  Auhl …")`。允许打印与 `examples/pack_peo` `lattice` 同形的一行。禁止在星形路径提到 Auhl、`CbmcGrow`、`PEO_GROW_TOL`、0.6 Å。

模块文档字符串同步。删除 `pack_star` 对 `_grow_tol()` 的调用；若该辅助函数不再被本文件引用则删除它。`PEO_GROW_TOL` 仍留在 `test_auhl_two_stars_push_off_converges` 和 Rust `examples/pack_peo` `auhl`。

`pack_ring` 不改算法：两生长器均对环抛 `RingTemplate`；示例仍先演示 `CbmcGrow` 具名拒绝，再选刚体 `GencanPack`。

### 文档

`docs/python/api-reference.md` 的 LatticeGrow 段删除 “Linear sp³ heavy-atom backbones only in v1 — branched and non-tetrahedral templates are named rejections.” 换为：

> Trees (including branched) are accepted; a cycle raises named `RingTemplate` `ValueError`; non-tetrahedral templates stay named rejections.

禁止：C3、BFS、「stars use lattice / stars are packed with LatticeGrow」、v1、linear-only。

`docs/python/guide/growth.md` 的 “Branched trees and rings” 必须同时成立：

- 化学仍由 molrs + molpy 构建。
- 两生长器都消费键图；molpack 从不按分子推断算法。
- 树（线性或分支）对 `CbmcGrow` **和** `LatticeGrow` 都合法。
- 走一遍 `pack_peo_topo.pack_star`：这是**一次显式 pick**（`LatticeGrow` @ 2.0 Å → 调用方 `GencanPack.with_restart` @ 2.0 Å）。
- `CbmcGrow` 仍是对等树生长器；Auhl 减核留在 `CbmcGrow` 路径。
- 环：两生长器都抛 `RingTemplate`；`pack_ring` 再选刚体 `GencanPack`。
- 禁止写 “stars are packed with LatticeGrow”。

### 测试

- **删除** `TestNamedErrors.test_lattice_refuses_star`。
- **新增** `TestNamedErrors.test_lattice_refuses_ring`：`make_ring` + `LatticeGrow.run`，`pytest.raises(ValueError, match="ring")`。
- **新增** `TestTinyPack.test_lattice_grows_star`：`make_star(2)` + `LatticeGrow` @ 2.0 Å、density 0.2；不 raise；`grown.natoms == star.n_atoms`；**无** push-off；不断言 `fdist == 0`。
- **保留** `test_cbmc_grows_one_star` 与 `test_auhl_two_stars_push_off_converges` 作为 P8 见证（后者继续硬编码 0.6 Å）。

### Reuse decision

- reuse `LatticeGrow` — 星形示例与 `test_lattice_grows_star` 的生长入口。
- reuse `GencanPack.with_restart` — 调用方 push-off；容差 2.0 Å。
- reuse `Target` / `PEO_H_RADIUS` / `TorsionPrior.three_state_from_c_inf`
- reuse `make_star` / `make_ring` / `pack_ring`
- reuse `CbmcGrow` — 对等树生长器；**不是** `pack_star` 配方。
- reuse `GrowError::RingTemplate`（Python `ValueError`，消息含 `ring`）
- reuse `examples/pack_peo` `lattice` 动词、`growth.md` 熔体格相样例、`test_lattice_then_seeded_push_off` — 2.0 Å 然后 `with_restart` 的同形。
- generalize `pack_star` — 探测拒绝 → 生长路径。
- generalize `test_lattice_refuses_star` — 改为 `test_lattice_grows_star`。
- new — 无。禁止 `StarGrow`、`pack_star_lattice`、第二套容差语义、示例内静默 fallback。

## Files to create or modify

- `python/tests/test_pack_peo_topo.py`
- `python/examples/pack_peo_topo.py`
- `docs/python/api-reference.md`
- `docs/python/guide/growth.md`
- ~~`regressions/lattice-branch-saw-02-docs.md` (new)~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）

## Tasks

- [x] Write failing unit tests for lattice star grow and ring refuse (`python/tests/test_pack_peo_topo.py` → `TestTinyPack.test_lattice_grows_star`, `TestNamedErrors.test_lattice_refuses_ring`); delete `test_lattice_refuses_star`; keep `test_cbmc_grows_one_star` and `test_auhl_two_stars_push_off_converges`
- [x] Implement `pack_star` in `python/examples/pack_peo_topo.py` as `LatticeGrow` @ 2.0 Å then `GencanPack.with_restart` @ 2.0 Å; delete lattice try/except and the `CbmcGrow` body; do not print Auhl; drop unused `_grow_tol`; keep `PEO_H_RADIUS` on `_target`
- [x] Replace the LatticeGrow v1 sentence in `docs/python/api-reference.md` with the trees / `RingTemplate` / non-tetrahedral acceptance set (no C3, no BFS, no “stars use lattice”)
- [x] Rewrite `docs/python/guide/growth.md` section “Branched trees and rings” as an explicit `pack_star` pick (`LatticeGrow` then caller-side `with_restart`); keep `CbmcGrow` a peer tree grower; rings raise `RingTemplate` on both growers and `pack_ring` picks rigid `GencanPack`
- [x] ~~Add regression example `regressions/lattice-branch-saw-02-docs.md` (public API only; hard-coded goldens, no third-party runtime)~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）
- [x] Verify owned files against the regression literals and the 2.0 Å / `RingTemplate` cases
- [x] Run full check + test suite

## Testing strategy

Python 单测所有权在 `python/tests/test_pack_peo_topo.py`。单文件绿：`tox -c python -e py -- tests/test_pack_peo_topo.py`（D-04：`uv run --directory python --group dev tox -e py` 可能因 molrs pin 无法解析）。格相星形单测只打 `LatticeGrow.run`，不跑 `pack_star` 全文。

- Happy：`make_star(2, seed=42)`，1 拷贝，`LatticeGrow` @ 2.0 Å、density 0.2；`grown.natoms == star.n_atoms`（`arm_length=2` 钉为 77）。无 push-off。
- Edge：`make_ring(4)` + `LatticeGrow.run` → `ValueError`，`match="ring"`。
- Peer：`test_cbmc_grows_one_star`；`test_auhl_two_stars_push_off_converges`（0.6 Å 然后 2.0 Å）。
- Domain：容差字面量 2.0 Å，禁止 0.6 Å 出现在格相星形测试或 `pack_star`。不断言格相 `fdist == 0`。
- ~~Regression：`regressions/lattice-branch-saw-02-docs.md` 钉字面量（日期 2026-09-05）。~~（`regressions/` 已于 2026-09-20 删除，无替代：golden 钉值不再是测试形式）

## Out of scope

- 任何 Rust / `python/src/` PyO3 / `src/grow/lattice/**` 改动（01-walk）。
- 公开新类型：`StarGrow`、`pack_star_lattice`。
- 文档或测试泄露 C3、BFS、子树 recoil。
- 改 `test_auhl_two_stars_push_off_converges` 或 Rust `examples/pack_peo` `auhl`。
- 改 `pack_ring` 的刚体 `GencanPack` 选择；不做格上闭环。
- 把 Auhl 0.6 Å 第一阶段抄到 `LatticeGrow`。
- 用 `Pipeline` 包装星形示例。
- 改 `docs/examples.md` / `docs/python/examples.md` 把星形绑定到 `LatticeGrow`。
- 力场、化学感知、`.inp` 新关键字。
