---
title: stl-region-02-bind — StlRegion 的 Python 绑定与文档
status: approved
created: 2026-09-05
chain: stl-region (02 of 02; predecessor stl-region-01-region)
---

# stl-region-02-bind — StlRegion 的 Python 绑定与文档

## Summary

01 落地后，Python 公开 `StlRegion.from_file(path, scale=1.0)`。`Target.with_restraint(stl)` 走逐原子通道，提取时 `RegionRestraint(inner)`。文档按**三个入口**写既有合同，不把 STL 写成一种 Force：GenCanPack 软惩罚、CbmcGrow propose 硬拒绝、LatticeGrow 屏蔽网格外格点。本段不改 `src/region/stl.rs`。

## Domain basis

不引入新方程。\(\varphi\) 与水密由 01 拥有。抬升 \(f=s_2\max(0,\varphi)^2\)。单位 Å。中心约束。

三个入口对**任何**几何 `AtomRestraint`（含 `InsideBoxRestraint` 与抬升后的 `StlRegion`）相同，不是 STL 特例：

- `GenCanPack`：二次外罚 + 梯度（软，中心可略越界）。
- `CbmcGrow`：`RestraintTable::violated` 在 propose 硬拒绝；`force_place` 不读 restraint，原子可在网格外并计 `softened`。
- `LatticeGrow`：Region ∩ lattice — 网格外的金刚石格点屏蔽；空交集 `GrowError::LatticeRegionEmpty`。

## Design

`python/src/region.rs`：`#[pyclass(name = "StlRegion")]` 包 `molpack::StlRegion`。只有 `from_file`。`from_file` 调用 01 的 `StlRegion::from_file`（`std::fs`，不启用 wheel 的 `io` feature）。禁止 Python 解析 STL。

`try_atom_builtin` 与 `extract_restraint` 各加一臂：`SharedAtomRestraint(Arc::new(RegionRestraint(inner)))`。`extract_collective_restraint` 不加臂。`StlRegion` **不**并入 `BuiltinRestraint` 的 `*Restraint` 名单；stub 里单独列出，作为 Region 走既有原子路径。

文档三页合同（guide/restraints.md、api-reference.md，必要时 growth.md 一句交叉）：上述三入口分裂；不要写「StlRegion/Force 只属于 GenCanPack」；不要从网格推断算法。

### Reuse decision

- generalize `try_atom_builtin` / `extract_restraint`
- reuse `RegionRestraint`, `StlRegion`, `Target.with_restraint`
- new `python/src/region.rs`
- new — 禁止 `StlRestraint`

## Files to create or modify

- `python/src/region.rs` (new)
- `python/src/constraint.rs`
- `python/src/lib.rs`
- `python/src/target.rs`
- `python/python/molpack/__init__.py`
- `python/python/molpack/molpack.pyi`
- `python/tests/test_region.py` (new)
- `docs/python/guide/restraints.md`
- `docs/python/api-reference.md`
- `regressions/stl-region-02-bind.md` (new)

## Tasks

- [ ] Write failing unit tests for StlRegion (`python/tests/test_region.py` → TestStlRegion)
- [ ] Implement PyStlRegion::from_file in `python/src/region.rs`
- [ ] Generalize try_atom_builtin and extract_restraint in `python/src/constraint.rs`
- [ ] Register StlRegion on the wheel (`lib.rs`, `__init__.py`, `.pyi`, target rustdoc)
- [ ] Document the three-entry split (Å, centre-only, GenCanPack soft / CbmcGrow propose / LatticeGrow Region ∩ lattice)
- [ ] Add regression example `regressions/stl-region-02-bind.md`
- [ ] Verify attach without TypeError and named parse errors
- [ ] Run full check + test suite

## Testing strategy

`python/tests/test_region.py`：tmp_path ASCII 水密单位立方体。`from_file` 成功；缺文件 OSError；非水密 ValueError。`Target.with_restraint` 与 `with_atom_restraint` 不 TypeError。不 `run()` 任何入口。`tox -c python -e py -- tests/test_region.py`。

## Out of scope

- src/region/stl.rs、grow、gencan
- StlRestraint、Python contains/f/fg、集体臂、`.inp` 关键字
