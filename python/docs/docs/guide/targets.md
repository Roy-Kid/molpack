# Targets

A `Target` describes one type of molecule to pack: its template
geometry, element symbols, and the number of copies to produce.
VdW radii are resolved automatically from element symbols via the
Bondi (1964) table.

## Construction

```python
from molpack import Target

target = Target(frame, count)
```

- `frame` — any object with `frame["atoms"]` returning a block.
  Supported sources:

  | Source | Element column | Block access |
  |--------|---------------|--------------|
  | `molrs.read_pdb(path)` | `"symbol"` | `.view()` |
  | `molrs.read_xyz(path)` | `"element"` | `.view()` |
  | `molpy.Frame` | `"element"` | `[]` |
  | plain dict | `"element"` | `[]` |

- `count` — number of copies to produce.

A display label is optional:

```python
target = Target(frame, count).with_name("water")
```

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
water = Target(frame, count=100).with_name("water")
```

## Read-only properties

```python
target.name        # Optional[str]
target.natoms      # number of template atoms
target.count       # requested copies
target.elements    # list[str]
target.radii       # list[float]
target.is_fixed    # True if placement is frozen (see below)
```

All builder methods are **immutable** — they return a new `Target`.

## Centering

The default is [`CenteringMode.AUTO`][molpack.CenteringMode] — free
targets are centered on their geometric center before packing; fixed
targets are kept in place. Override explicitly:

```python
from molpack import CenteringMode

target = target.with_centering(CenteringMode.CENTER)  # always center
target = target.with_centering(CenteringMode.OFF)     # keep input coords
```

## Fixed placement

Pin a target at a specific location (e.g. a reference protein):

```python
from molpack import Angle

target = target.fixed_at([10.0, 20.0, 30.0])

# optional Euler orientation — three Angle values in Packmol's
# eulerfixed convention
target = (
    target.fixed_at([10.0, 20.0, 30.0])
    .with_orientation((
        Angle.from_degrees(0.0),
        Angle.from_radians(1.57),
        Angle.ZERO,
    ))
)
```

Fixed targets are excluded from the optimizer but still contribute to
distance exclusion against other species.

## Rotation bounds

Restrict the rotational search window about each axis:

```python
from molpack import Angle, Axis

target = (
    target
    .with_rotation_bound(Axis.X, Angle.from_degrees(0.0),  Angle.from_degrees(15.0))
    .with_rotation_bound(Axis.Y, Angle.from_degrees(90.0), Angle.from_degrees(10.0))
    .with_rotation_bound(Axis.Z, Angle.from_degrees(0.0),  Angle.from_degrees(5.0))
)
```

`Angle` makes units explicit: use `Angle.from_degrees(...)` or
`Angle.from_radians(...)` — raw floats are rejected.

## Attaching restraints

### All atoms of the target

```python
from molpack import InsideBoxRestraint

target = target.with_restraint(
    InsideBoxRestraint([0, 0, 0], [40, 40, 40])
)
```

Stack multiple restraints by calling `.with_restraint()` again:

```python
target = (
    target
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
    .with_restraint(OutsideSphereRestraint([20, 20, 20], 5.0))
)
```

### A subset of atoms

```python
from molpack import AbovePlaneRestraint, BelowPlaneRestraint

target = target.with_atom_restraint(
    [30, 31],                                     # 0-based Rust-native indices
    BelowPlaneRestraint([0.0, 0.0, 1.0], 2.0),
)
```

!!! note "0-based indexing"
    `with_atom_restraint` uses **0-based** indices, matching Rust
    convention. If you are porting from a Packmol `.inp` file (which
    uses 1-based indices), subtract 1 at the call site.

## Per-target solver budget

Override the maximum perturbation budget for this target:

```python
target = target.with_perturb_budget(50)  # default: derived from count
```

Useful when one species is significantly harder to place than the rest.
