---
title: stl-region-01-region — 封闭三角网格作为 Region
status: approved
created: 2026-09-05
chain: stl-region (01 of 02; successor stl-region-02-bind)
---

# stl-region-01-region — 封闭三角网格作为 Region

## Summary

调用方把 watertight 三角网格变成普通 `Region`：`StlRegion` 经 `from_file` → `from_bytes` →（缩放成 Å）→ `from_triangles` 构造，实现与盒/球相同的 `contains` / `signed_distance` / `signed_distance_grad` / `bounding_box`，再经既有 `RegionRestraint` 接到 `Target.with_restraint`。不新增 `StlRestraint`。表面（无符号距离 `< 1e-9` Å）视为在内且 `φ = 0`；否则半开 Möller–Trumbore 偶奇 + Eberly 最近点，符号负内正外。非法网格在构造期具名拒绝。本段只改 `src/region/` 与 crate 根 re-export。

## Domain basis

`contains(x) ⇔ signed_distance(x) ≤ 0`。单位 Å。只测原子中心。网格是有限 \(\mathbb{R}^3\) 中的闭 2-流形，查询点不按 PBC 折回。

无符号距离 = Eberly 最近点到任一三角形（七区 Voronoi：面 / 三边 / 三顶点）。doi 文献：Eberly, *Distance Between Point and Triangle in 3D*, Geometric Tools (1999), https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

\(\varepsilon = 10^{-9}\) Å：\(d<\varepsilon \Rightarrow \varphi=0\)（在内）。否则 \(\varphi=\sigma d\)，\(\sigma\) 来自固定无理射线 \(\hat r=(1,\sqrt{2},\pi)/\|\cdot\|\) 的偶奇命中。求交：Möller & Trumbore, *J. Graphics Tools* **2**, 21–28 (1997), doi:10.1080/10867651.1997.10487468；半开重心 \(u\ge 0,v\ge 0,u+v<1\)，\(t>\varepsilon\)。嵌套闭壳按偶奇成洞，不用绕数/伪法向。

梯度：\(d<\varepsilon\) 为零；否则 \(\sigma(x-c)/d\)（指向外侧）。`RegionRestraint` 只在 \(\varphi>0\) 使用。

教科书立方体 \([0,1]^3\) Å、12 三角、`scale=1`：中心 \(\varphi=-0.5\)；面 \((1,0.5,0.5)\) \(\varphi=0\) 且 contains；\((2,0.5,0.5)\) \(\varphi=+1\)；\((3,3,3)\) \(\varphi=2\sqrt{3}\)。

## Design

机械迁徙 `src/region.rs` → `src/region/mod.rs`（公开符号不变）。`mod stl` **私有**；`pub use stl::{StlRegion, StlError}`。`stl.rs` 不 import `restraint`。

构造栈（一个身份）：`from_file(path, scale) → from_bytes(bytes, scale) → 顶点×scale → from_triangles(Å)`。三者都返回 `Result<Self, StlError>`，跑同一套门。`from_triangles` **没有** scale（已是 Å）。加载器不保留路径、solid 名、STL 法向、attribute。

顶点焊接（构造内部，不让调用方按约定预焊）：缩放后按 `f64::to_bits` 三元组精确合并；然后每条无向边必须恰好两条反向半边。仍开/非流形/退化/空/非有限 → `StlError`。不按容差 snap（政策写死：bitwise 重合才是同一顶点）。

存储：私有 `Vec`（或 boxed slice）。禁止公开 `Arc`、三角数组、半边、射线方向。`Clone` 可以拷贝三角；热路径共享走既有 `Arc<dyn AtomRestraint>`。

`StlError`：`Io` / `Parse` / `Empty` / `InvalidScale` / `DegenerateTriangle { index }` / `NotWatertight { unpaired }`。`Clone`。不进 `PackError`。不挂 `io` feature。

`declared_cell` 恒 `None`。`bounding_box` 为顶点 AABB（Å）。

CBMC Force / LatticeGrow：本段不改 `grow/`；不声称所有求解器构造性限制在网格内。

### Reuse decision

- reuse `Region`, `RegionRestraint`, `Aabb`, `RegionExt::into_restraint`
- new `StlRegion` / `StlError` / STL 解析（无既有网格谓词）
- new — 禁止 `StlRestraint`

## Files to create or modify

- `src/region.rs` (delete; content moves)
- `src/region/mod.rs` (new)
- `src/region/stl.rs` (new)
- `src/lib.rs`
- `regressions/stl-region-01-region.md` (new)

## Tasks

- [x] Move `src/region.rs` to `src/region/mod.rs` (public symbols unchanged)
- [ ] Write failing unit tests for `StlRegion` / `StlError` (`src/region/stl.rs` `#[cfg(test)]`)
- [ ] Implement `StlRegion`, `StlError`, and STL parse in `src/region/stl.rs`; private `mod stl` from `mod.rs`
- [ ] Add rustdoc per rustdoc with units (Å, scale, φ sign, ε = 1e-9 Å)
- [ ] Re-export `StlRegion` and `StlError` from `src/lib.rs` (crate root, prelude, table)
- [ ] Add regression example `regressions/stl-region-01-region.md`
- [ ] Verify unit-cube SDF, nested even-odd cavity, named rejects
- [ ] Run full check + test suite

## Testing strategy

In-module `src/region/stl.rs`。`cargo test -p molcrafts-molpack --lib --tests -- region::stl`。迁徙：`region::` 旧测试仍绿。立方体金标见 Domain basis。ASCII 与 length-matched binary 与 `from_triangles` 一致。`RegionRestraint` 中心 `f==0`、面外 `f==scale2`。

## Out of scope

- python/、docs/python/、`.inp`（02-bind）
- grow/、gencan/、BVH、容差焊接、绕数、PBC 折网格、半径膨胀
- `StlRestraint`、`PackError` STL 变体、新 crate 依赖
