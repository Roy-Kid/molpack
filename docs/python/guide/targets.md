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

- `frame` — a `molrs.core.Frame` (`molpy.Frame` is the same class), resolved
  zero-copy via its FFI capsule. Element symbols come from the `"element"`
  atom column, which every molrs reader (`molrs.io.read_pdb`,
  `molrs.io.read_xyz`, …) writes; a frame built in memory must carry it too.

- `count` — number of copies to produce.

A display label is optional:

```python
target = Target(frame, count).with_name("water")
```

Build a frame in memory (no PDB file) with `molrs.core.Frame`:

```python
import molrs
import numpy as np

frame = molrs.core.Frame({
    "atoms": {
        "x": np.array([0.00,  0.96, -0.24]),
        "y": np.array([0.00,  0.00,  0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
})
water = Target(frame, count=100).with_name("water")
```

## Read-only properties

```python
target.name           # Optional[str]
target.natoms         # number of template atoms
target.count          # requested copies
target.elements       # list[str]
target.radii          # list[float]
target.special_bonds  # list[float] — intramolecular skip table
target.is_fixed       # True if placement is frozen (see below)
```

All builder methods are **immutable** — they return a new `Target`.

## Centering

The default is [`CenteringMode.AUTO`](../api-reference.md#centeringmode) — free
targets are centered on their geometric center before packing; fixed
targets are kept in place. Override explicitly:

```python
from molpack import CenteringMode

target = target.with_centering(CenteringMode.CENTER)  # always center
target = target.with_centering(CenteringMode.OFF)     # keep input coords
```

## Fixed placement

Pin a target at a specific location (e.g. a fixed reference molecule):

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
import molrs

target = target.with_restraint(
    molrs.core.Cuboid([0, 0, 0], [40, 40, 40])   # a molrs region: origin, lengths
)
```

Stack multiple restraints by calling `.with_restraint()` again:

```python
target = (
    target
    .with_restraint(molrs.core.Cuboid([0, 0, 0], [40, 40, 40]))
    .with_restraint(~molrs.core.Sphere([20, 20, 20], 5.0))
)
```

### A subset of atoms

```python
target = target.with_atom_restraint(
    [30, 31],                                     # 0-based Rust-native indices
    molrs.core.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 2.0]),   # z <= 2
)
```

!!! note "0-based indexing"
    `with_atom_restraint` uses **0-based** indices, matching Rust
    convention. If you are porting from a Packmol `.inp` file (which
    uses 1-based indices), subtract 1 at the call site.

## Packing radii

The packer separates two atoms by the sum of their packing radii;
without an override every atom uses the global `tolerance / 2`. Van der
Waals radii from the source file are not used as packing radii.

```python
target = target.with_radius(2.0)                    # every atom
target = target.with_atom_radius(h_indices, 0.85)   # then explicit hydrogen
```

All-atom chains with explicit hydrogen keep the default depth-3 skip
table and shrink hydrogen here (~0.85 Å). Indices are **0-based**.

## Intramolecular skip table

Atom pairs close along the chain are exempt from the hard core. The
table is per-target data, not an engine knob:

```python
target.special_bonds  # [0.0, 0.0, 0.0, 1.0] by default (depth 3)
cg = target.with_special_bonds([0.0, 0.0, 1.0])  # CG depth 2
```

Empty, non-finite, or out-of-range weights raise `ValueError`.
Fractional weights are stored and refused later at `CbmcGrow.run`. See
[Chain growth](growth.md) for the all-atom vs CG recipe.

## Per-target solver budget

Override the maximum perturbation budget for this target:

```python
target = target.with_perturb_budget(50)  # default: derived from count
```

Useful when one species is significantly harder to place than the rest.
