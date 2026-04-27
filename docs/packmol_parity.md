# Packmol parity

molpack tracks the original Packmol's behavior (`L. Martínez et al.,
J. Comput. Chem. 2009`). Coordinates are not byte-identical with the
same seed, but functional equivalence is enforced on five canonical
workloads.

## What is matched

**Objective structure**

- Geometric restraint penalties — equivalent of `comprest` / `gwalls`.
- Minimum-distance overlap penalty — equivalent of `computef` /
  `fparc`.

**Optimization workflow**

- Initialization with constraint-only fitting (`initial` / `restmol` /
  `swaptype`).
- Phased main optimization — per-type pre-compaction, then all-types.
- `movebad` heuristic for stalled molecules.
- Radius-scaling (`radscale`) decay across each phase.
- GENCAN / SPG / CG inner solver chain.
- Precision gate on `fdist` (overlap) and `frest` (restraint
  violation).

**Restraint vocabulary used by the canonical examples**

- `inside box`, `inside sphere`, `outside sphere`
- `above plane` / `below plane`
- fixed molecule placement

**Periodicity**

- `pbc X Y Z` and `pbc X0 Y0 Z0  X1 Y1 Z1`, mirroring Packmol's
  `getinp.f90`. The packer's cell grid is built directly from the
  declared box.

**Script parser strictness**

- Unknown top-level keywords are rejected via
  `ScriptError::UnknownKeyword`. A silently dropped `pbc` previously
  triggered a 42 GB cell-grid allocation — strict parsing prevents
  that class of failure.

**Determinism**

- Explicit seeds; identical seed values are used for paired Packmol
  vs molpack runs in the regression suite.

## Verification

### Batch example validation

`tests/examples_batch.rs` runs all five canonical workloads. Marked
`#[ignore]` because the run is expensive — invoke explicitly:

```bash
cargo test -p molcrafts-molpack --release --test examples_batch -- --ignored
```

The test asserts:

- atom and molecule counts match the expanded target specs;
- XYZ output is structurally sound;
- quantified violation metrics stay within tolerance / precision.

### Violation metrics

Each side reports `max_distance_violation`, `max_constraint_penalty`,
`violating_pairs`, and `violating_atoms`. Both packers must satisfy the
same thresholds.

## Accepted differences

- Packed coordinates are not bit-identical with Packmol even at the
  same seed.
- Acceptance is functional: same restraints, same conflict criteria,
  same magnitude of tolerance / precision, comparable violation metrics.
