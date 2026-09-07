# Regression — stl-region-02-bind

Public API. 2026-09-05. No third-party runtime.

- `StlRegion.from_file(path, scale=1.0)` — Å per file unit.
- `Target.with_restraint(stl)` lifts per-atom via `RegionRestraint` (not collective).
- No `StlRestraint`. No Python `f`/`fg` on `StlRegion`.
- Three-entry split: GenCanPack soft; CbmcGrow propose hard / Force may be outside; LatticeGrow masks diamond sites outside the region (Region ∩ lattice). Empty intersection is `GrowError::LatticeRegionEmpty`.
