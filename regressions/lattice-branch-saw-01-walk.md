# Regression — lattice-branch-saw-01-walk

Public API only. In-crate fixture, 2026-09-05. No third-party runtime.

```rust
use molpack::{LatticeGrow, Target, TorsionPrior};

// Tetrahedral 5-atom star: centre at origin, leaves along diamond T,
// bond 1.53 Å. Two copies, seed 7, cubic box 20 Å.
let grown = LatticeGrow::new(TorsionPrior::Uniform)
    .with_seed(7)
    .with_tolerance(2.0)
    .with_periodic_box([0.0; 3], [20.0; 3], [true; 3])
    .run(&[Target::new(star_frame, 2)], 60)
    .unwrap();
assert_eq!(grown.natoms(), 10);
// centre–leaf bonds 1.53 ± 1e-6 (template InternalTree geometry)
```

Owning tests in `tests/grow.rs`:

- `lattice_grow_tetrahedral_star_completes` — `Ok`, `natoms == 10`
- `lattice_grow_tetrahedral_comb_completes` — `Ok`, `natoms == 24`
- `lattice_grow_rejects_degree_gt_4` — `NonTetrahedralTemplate`, no “branched staged”
- `lattice_grow_bead_chain_constructive` — linear `k = 1`: `fdist` bitwise 0, `softened == 0`, bonds 1.53 Å, same-seed bitwise
