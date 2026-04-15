# Targets

A `Target` describes one type of molecule to pack: its template
geometry, element symbols, and the number of copies to produce.
VdW radii are resolved automatically from element symbols via the
Bondi (1964) table.

## Construction

```python
from molpack import Target

target = Target(name, frame, count)
```

- `name`  — label used in diagnostics and output.
- `frame` — any object with `frame["atoms"]` returning a block.
  Supported sources:

  | Source | Element column | Block access |
  |--------|---------------|--------------|
  | `molrs.read_pdb(path)` | `"symbol"` | `.view()` |
  | `molrs.read_xyz(path)` | `"element"` | `.view()` |
  | `molpy.Frame` | `"element"` | `[]` |
  | plain dict | `"element"` | `[]` |

- `count` — number of copies to produce.

Plain dict example (no I/O library needed):

```python
import numpy as np

frame = {
    "atoms": {
        "x": np.array([0.00,  0.96, -0.24]),
        "y": np.array([0.00,  0.00,  0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
}
water = Target("water", frame, count=100)
```

## Read-only properties

```python
target.natoms      # number of template atoms
target.count       # requested copies
target.elements    # list[str]
target.is_fixed    # True if placement is frozen (see below)
```

All builder methods are **immutable** — they return a new `Target`.

## Renaming

```python
target = target.with_name("water")   # update label after construction
```

## Centering

By default the template is re-centered on its mean position before
packing (Packmol convention). Control explicitly:

```python
target = target.with_center()        # force centering on
target = target.without_centering()  # keep template origin
```

## Fixed placement

Pin a target at a specific location (e.g. a reference protein):

```python
target = target.fixed_at([10.0, 20.0, 30.0])
target = target.fixed_at_with_euler(
    position=[10.0, 20.0, 30.0],
    euler=[0.0, 1.57, 0.0],          # radians
)
```

Fixed targets are excluded from the optimizer but still contribute to
distance exclusion against other species.

## Rotation bounds

Restrict the rotational search window about each axis (degrees):

```python
target = target.constrain_rotation_x(center_deg=0.0,  half_width_deg=15.0)
target = target.constrain_rotation_y(center_deg=90.0, half_width_deg=10.0)
target = target.constrain_rotation_z(center_deg=0.0,  half_width_deg=5.0)
```

## Attaching restraints

### All atoms of the target

```python
from molpack import InsideBox

target = target.with_restraint(InsideBox([0, 0, 0], [40, 40, 40]))
```

Stack multiple restraints by calling `.with_restraint()` again:

```python
target = (
    target
    .with_restraint(InsideBox([0, 0, 0], [40, 40, 40]))
    .with_restraint(OutsideSphere(5.0, [20, 20, 20]))
)
```

### A subset of atoms

```python
from molpack import AbovePlane, BelowPlane

target = target.with_restraint_for_atoms(
    [31, 32],                          # 1-based Packmol atom indices
    BelowPlane([0.0, 0.0, 1.0], 2.0),
)
```

!!! note "1-based indexing"
    `with_restraint_for_atoms` uses Packmol's 1-based atom indices.
    Passing `0` or an index `> natoms` raises `ValueError`.

## Per-target solver budget

Override the maximum optimizer iterations for this target:

```python
target = target.with_maxmove(50)  # default: derived from molecule size
```

Useful when one species is significantly harder to place than the rest.
