# API Reference

Import surface:

```python
from molpack import (
    Target, Molpack, PackResult,
    InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane,
)
```

---

## `Target`

Molecule-type specification: template geometry, element symbols, copy
count. VdW radii are looked up automatically from element symbols
(Bondi 1964). Immutable — builder methods return new instances.

**Constructor**

```python
Target(name: str, frame, count: int)
```

- `name`  — label used in diagnostics and output.
- `frame` — any object with `frame["atoms"]` returning a block with
  columns `"x"`, `"y"`, `"z"`, and `"element"` (or `"symbol"` for
  `molrs` PDB frames). Accepts `molrs.Frame`, `molpy.Frame`, or a
  plain dict.
- `count` — number of copies to produce.

**Metadata builders**

- `.with_name(name: str) -> Target`
- `.with_center() -> Target`
- `.without_centering() -> Target`

**Placement**

- `.fixed_at(position: [x, y, z]) -> Target`
- `.fixed_at_with_euler(position, euler) -> Target` — `euler` in radians.

**Rotation bounds** (degrees)

- `.constrain_rotation_x(center_deg, half_width_deg) -> Target`
- `.constrain_rotation_y(center_deg, half_width_deg) -> Target`
- `.constrain_rotation_z(center_deg, half_width_deg) -> Target`

**Restraints**

- `.with_restraint(restraint) -> Target` — attach to all atoms; call
  multiple times to stack restraints.
- `.with_restraint_for_atoms(indices: list[int], restraint) -> Target`
  — 1-based Packmol atom indices.

**Solver budget**

- `.with_maxmove(maxmove: int) -> Target`

**Properties (read-only)**

- `.natoms : int`
- `.count : int`
- `.elements : list[str]`
- `.is_fixed : bool`

---

## `Molpack`

Orchestrator for the three-phase GENCAN optimizer.

**Constructor**

```python
Molpack(
    tolerance=2.0,
    precision=0.01,
    maxit=20,
    nloop0=0,
    sidemax=1000.0,
    movefrac=0.05,
    movebadrandom=False,
    disable_movebad=False,
    progress=True,
)
```

**Builder methods** (all immutable)

- `.with_tolerance(t: float)`
- `.with_precision(p: float)`
- `.with_maxit(n: int)`
- `.with_nloop0(n: int)`
- `.with_sidemax(s: float)`
- `.with_movefrac(f: float)`
- `.with_movebadrandom(enabled: bool)`
- `.with_disable_movebad(disabled: bool)`
- `.with_pbc(min: [x,y,z], max: [x,y,z])`
- `.with_pbc_box(lengths: [x,y,z])`
- `.with_progress(enabled: bool)`

**Run**

```python
.pack(targets: list[Target], max_loops: int = 200, seed: int | None = None) -> PackResult
```

---

## `PackResult`

Read-only output container.

**Properties**

- `.positions : ndarray (N, 3) float64` — packed coordinates.
- `.elements  : list[str]`
- `.natoms    : int`
- `.converged : bool`
- `.fdist     : float`
- `.frest     : float`

---

## Restraints

All restraint classes are immutable.

### `InsideBox(min, max)`

Axis-aligned box. `min`, `max` are 3-element sequences.

### `InsideSphere(radius, center)`

Closed ball.

### `OutsideSphere(radius, center)`

Complement of closed ball.

### `AbovePlane(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \ge d\}$.

### `BelowPlane(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \le d\}$.

---

## Exceptions

`molpack` raises standard Python exceptions:

- `ValueError` — invalid atom index (must be 1-based and in $[1, N]$),
  or mismatched array lengths when building a frame manually.
- `RuntimeError` — optimizer-level failures (zero targets, invalid PBC
  box, etc.). Message is the Rust `PackError::Display` rendering.
- `TypeError` — wrong argument type passed to restraint builders.
