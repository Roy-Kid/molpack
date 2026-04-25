# Python Targets

`Target` is the heart of the Python API. It describes one molecule type: the
template coordinates, element symbols, and the number of copies the solver
should place.

## Start With the Simple Mental Model

Think of a `Target` as a species description, not as one already-placed
molecule. `count=100` means "make one hundred copies of this template while
solving", not "here are one hundred Python objects".

```python
from molpack import Target

target = Target(frame, count=100)
```

Supported frame-like inputs include `molrs.Frame`, `molpy.Frame`, and ordinary
Python dicts with an `atoms` block containing `x`, `y`, `z`, and `element`
columns.

## Give It a Name

Naming a target is optional, but it pays off fast when you inspect logs or
debug mixed systems.

```python
water = Target(frame, count=100).with_name("water")
```

All builder methods are immutable. They return a new `Target` instead of
modifying the old one in place.

## Attach Geometry

Every mobile target needs at least one geometric restraint.

```python
from molpack import InsideBoxRestraint

water = water.with_restraint(
    InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])
)
```

You can stack restraints with repeated `.with_restraint(...)` calls or attach a
restraint to only some atoms with `.with_atom_restraint(indices, restraint)`.
Atom indices are `0`-based.

## Centering, Fixed Placement, and Rotation Bounds

Free targets are usually centered automatically before packing. Fixed targets
are the special case: they stay where you pin them.

```python
from molpack import Angle, Axis, CenteringMode

target = target.with_centering(CenteringMode.AUTO)
target = target.fixed_at([10.0, 20.0, 30.0])
target = target.with_orientation((
    Angle.ZERO,
    Angle.from_degrees(90.0),
    Angle.ZERO,
))
target = target.with_rotation_bound(
    Axis.Z,
    Angle.from_degrees(0.0),
    Angle.from_degrees(5.0),
)
```

Use fixed placement when one structure should act as an immobile reference, for
example a protein in a solvation setup.

## A Plain-Dict Example

You do not need an I/O package for the smallest examples:

```python
import numpy as np
from molpack import InsideBoxRestraint, Target

frame = {
    "atoms": {
        "x": np.array([0.00, 0.96, -0.24]),
        "y": np.array([0.00, 0.00, 0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
}

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0.0, 0.0, 0.0], [30.0, 30.0, 30.0]))
)
```

## Read-Only Properties

These are the properties you will reach for most often:

- `target.name`
- `target.natoms`
- `target.count`
- `target.elements`
- `target.radii`
- `target.is_fixed`

## Takeaway

If you remember one thing, make it this: `Target` is where molecule-specific
facts live. Copy count, fixed placement, rotation limits, and target-local
restraints all belong here, not on the packer.
