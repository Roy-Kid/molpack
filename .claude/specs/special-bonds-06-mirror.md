---
title: special-bonds-06-mirror — Python Target.with_special_bonds 与文档转向
status: done
created: 2026-09-04
depends_on: [special-bonds-05-target]
chain: special-bonds（06 of 6；Python 绑定 + 文档层）
grilled: true
---

# special-bonds-06-mirror — Python Target.with_special_bonds 与文档转向

## Summary

把 05 已经收在 Rust `Target` 上的特殊键表镜像到 Python wheel：`Target.with_special_bonds(list[float])` 是唯一入口，默认表仍是深度 3 的 `[0, 0, 0, 1]`。引擎旋钮由 05 解钩；本 spec 验收它已经不在，再补 `Target.with_special_bonds`。`PackResult` 带上 04 的嵌套残差 `intra: IntraResidual { scored, exempted }`。文档把「显式氢用 `Target.with_atom_radius` 给小半径」写成全原子链的一等建议（2026-09-04 实测：深度 3、H = 0.85 Å → 142 轮 / 0.3 s），`[0,0,0,0,0,1]` 只作为残差报告暴露的逃生阀。

## Domain basis

本 spec 不改生长几何。权重无量纲；残差与半径的单位是 Å。出处与 05 相同。实测表（molpack 真机，dp5 PEO，2026-09-04）硬编码进文档与 regression，测试期不重跑该熔体。

## Design

**05 公开面（逐字遵守）：**

```
impl Target {
    pub fn with_special_bonds(self, table: BondDistanceWeights) -> Self
}
pub special_bonds: BondDistanceWeights
// GrowError::NonBinarySpecialBond wrapped in PackError::Grow
// 05 already unhooks python/src/entry.rs, molpack.pyi, test_grow.py
```

**Python builder：** `&self` + clone `inner`，调用 05 的 `Self` 方法。`validate_special_bonds` 紧挨 `check_positive`，只做 `BondDistanceWeights::new` 的 marshalling：非空、有限、`0≤w≤1` → 否则 `ValueError`。**不**复制 `binary_violation`：分数 0.5 存进 `Target.special_bonds`，与 Rust `Self` 一致；非二值在 `CbmcGrow.run` / `LatticeGrow.run` 上以既有 `PackError::Grow`（`GrowError::NonBinarySpecialBond`）出现。`test_target.py` 不 `run()`；0.5 的具名错误放在 `test_grow.py` `TestNamedErrors`，紧挨 `test_grow_without_bonds`。getter 读 `inner.special_bonds.as_slice()`。不暴露 `BondDistanceWeights` pyclass。不改 `helpers.rs`。

**Files 不含 `python/src/entry.rs`。** 05 已经解钩。06 用 `hasattr` / ripgrep 验收缺席。

**残差：** `PyIntraResidual { scored, exempted }` 嵌套，`result.intra.scored` / `.exempted`。注册于 `lib.rs` 与 `__init__.py`。禁止 `min_intra_*`。不重算。

**stub：** 列出每一个已落地 Target builder，外加 `with_special_bonds` / `special_bonds`。

**测试所有权：** `test_target.py` 不 `run()`。`test_grow.py` 验收引擎方法不在。`test_pack_result.py` 转发，不从 positions 重推。

**文档：** `growth.md` AA/CG 整段改写为 `Target.with_special_bonds`；氢半径一等；`[0,0,0,0,0,1]` 只出现在读残差段。`docs/python/` 内 `with_exclusion_depth` 命中数为 0。

**Reuse decision**

- `PyTabulatedPlane` `Vec<NpF>` — **pattern**。
- Target builders / `check_positive` — **reuse**。
- entry.rs 五处引擎旋钮 — **already 05**。
- PackResult getters — **reuse**。
- `PyStageInfo` — **pattern** for nested IntraResidual。
- `helpers.rs` — **reuse as-is**。

## Files to create or modify

- `python/src/target.rs`
- `python/src/result.rs`
- `python/src/lib.rs`
- `python/python/molpack/molpack.pyi`
- `python/python/molpack/__init__.py`
- `python/tests/test_target.py` (new)
- `python/tests/test_grow.py`
- `python/tests/test_pack_result.py`
- `docs/python/guide/growth.md`
- `docs/python/api-reference.md`
- `docs/python/guide/targets.md`
- `docs/python/guide/packer.md`
- `docs/python/getting-started.md`
- `docs/concepts.md`
- `docs/extending.md`
- `docs/architecture.md`
- `regressions/special-bonds-06-mirror.md` (new)

## Tasks

- [x] Write failing unit tests for `Target.with_special_bonds` in `python/tests/test_target.py` (default table, immutability, CG table round-trip, empty/non-finite/out-of-range → ValueError, fractional 0.5 stores and getter round-trips; no run())
- [x] Implement `PyTarget.with_special_bonds` and `special_bonds` getter (`validate_special_bonds` sibling of `check_positive`, including binary 0/1)
- [x] Write failing tests in `test_pack_result.py` (nested IntraResidual forwarding; bonded diatomic for empty scored class `+∞`) and `test_grow.py` (`hasattr(..., "with_exclusion_depth") is False`; fractional 0.5 raises at `CbmcGrow.run` as GrowError::NonBinarySpecialBond)
- [x] Implement `PyIntraResidual` + `PackResult.intra`; register in `lib.rs`; export from `__init__.py`
- [x] Sync molpack.pyi (every live Target builder plus new table and IntraResidual)
- [x] Update docs (growth.md AA/CG rewrite; H-radius first-class; no with_exclusion_depth in docs/python/)
- [x] Add regression example `regressions/special-bonds-06-mirror.md`
- [x] Run `tox -c python -e py` plus python crate fmt/check and `cargo doc -p molcrafts-molpack --no-deps`

## Testing strategy

不重跑 peo-tg 熔体。默认表观测 05 公开字段。0.5 在 Target 上可存储，在 CbmcGrow.run 上具名失败。残差用有键双原子断言空计分类 +inf；`_make_tiny_pack` 只测名字/类型。不从 positions 重算 min。

## Out of scope

- 不改 Rust Target / 生长引擎（05）。不重删 entry.rs 引擎旋钮。
- 不重开 IntraResidual 嵌套命名。不发明 min_intra_* / resolved_special_bonds。
- 不改 helpers.rs / interop.rs。仓外脚本不在 Files。
- 不改 CLI / `.inp`。P10 不碰 capsule。
