# API Reference

Import surface:

```python
from molpack import (
    # Target + packer
    Target, Packer, PackResult,
    # Restraints
    InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane,
    # Composite
    MoleculeConstraint,
)
```

---

## `Target`

Molecule-type specification: template coords, radii, elements, copy
count. Immutable — builder methods return new instances.

**Constructor**

- `Target.from_coords(positions, radii, count, elements=None)` →
  `Target`
    - `positions : ndarray (N, 3) float64`
    - `radii     : ndarray (N,)   float64`
    - `count     : int`
    - `elements  : list[str] | None` — default `"X"` per atom.

**Metadata builders**

- `.with_name(name: str) -> Target`
- `.with_center() -> Target`
- `.without_centering() -> Target`

**Placement**

- `.fixed_at(position: [x, y, z]) -> Target`
- `.fixed_at_with_euler(position, euler) -> Target` — `euler` in radians.

**Rotation bounds**

- `.constrain_rotation_x(center_deg, half_width_deg) -> Target`
- `.constrain_rotation_y(center_deg, half_width_deg) -> Target`
- `.constrain_rotation_z(center_deg, half_width_deg) -> Target`

**Restraints**

- `.with_constraint(restraint) -> Target`
- `.with_constraint_for_atoms(indices: list[int], restraint) -> Target`
  — 1-based indexing.

**Solver budget**

- `.with_maxmove(maxmove: int) -> Target`

**Properties (read-only)**

- `.natoms : int`
- `.count : int`
- `.elements : list[str]`
- `.is_fixed : bool`

---

## `Packer`

Orchestrator for the three-phase GENCAN optimizer.

**Constructor**

- `Packer(tolerance=2.0, precision=0.01)`

**Builder methods** (all immutable)

- `.with_tolerance(t: float)`
- `.with_precision(p: float)`
- `.with_maxit(n: int)`
- `.with_nloop0(n: int)`
- `.with_sidemax(s: float)`
- `.with_movefrac(f: float)`
- `.with_movebadrandom(enabled: bool)`
- `.with_disable_movebad(disabled: bool)`
- `.with_pbc(min: [x, y, z], max: [x, y, z])`
- `.with_pbc_box(lengths: [x, y, z])`
- `.with_progress(enabled: bool)`

**Run**

- `.pack(targets: list[Target], max_loops: int, seed: int | None = None)
   -> PackResult`

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

All restraint classes are immutable and support `.and_(other)` to
compose into a `MoleculeConstraint`.

### `InsideBox(min, max)`

Axis-aligned box. `min`, `max` are 3-element sequences.

### `InsideSphere(radius, center)`

Closed ball.

### `OutsideSphere(radius, center)`

Complement of closed ball.

### `AbovePlane(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \ge d\}$.
`normal` is a 3-element sequence; the packer does not auto-normalise.

### `BelowPlane(normal, distance)`

Half-space $\{\mathbf{x} : \mathbf{n}\cdot\mathbf{x} \le d\}$.

### `MoleculeConstraint`

A bundle of restraints. Construct by calling `.and_()` on any
restraint; accepts nested `MoleculeConstraint` arguments. Pass to
`Target.with_constraint` like any single restraint.

---

## Exceptions

`molpack` raises standard Python exceptions from `pack()`:

- `ValueError` — shape mismatch, invalid atom index (must be 1-based
  and in $[1, N]$), invalid PBC box.
- `RuntimeError` — optimizer-level failures (produced from
  `PackError::Display`).
- `TypeError` — wrong argument type passed to restraint builders.

No custom exception hierarchy is exposed; consume the message string
for diagnostics.
