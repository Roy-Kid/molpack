# Regression — stl-region-01-region

Public API. 2026-09-05. No third-party runtime.

```rust
use molpack::{Region, RegionExt, StlRegion};

// 12-triangle unit cube [0,1]³ Å, scale already Å.
let cube = StlRegion::from_triangles(&unit_cube_tris).unwrap();
assert!(cube.contains(&[0.5, 0.5, 0.5]));
assert_eq!(cube.signed_distance(&[0.5, 0.5, 0.5]), -0.5);
assert_eq!(cube.signed_distance(&[2.0, 0.5, 0.5]), 1.0);
assert_eq!(cube.into_restraint().f(&[0.5, 0.5, 0.5], 1.0, 1.0), 0.0);
```

Owning tests: `src/region/stl.rs` `#[cfg(test)]`.
