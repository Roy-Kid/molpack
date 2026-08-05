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

- `frame` — a `molrs.Frame` or `molpy.Frame`, resolved zero-copy via its
  FFI capsule. Supported sources:

  | Source | Element column |
  |--------|---------------|
  | `molrs.io.read_pdb(path)` | `"symbol"` |
  | `molrs.io.read_xyz(path)` | `"element"` |
  | `molrs.Frame.from_dict({"blocks": {"atoms": {...}}, "meta": {}})` | `"element"` |
  | `molpy.Frame` | `"element"` |

- `count` — number of copies to produce.

A display label is optional:

```python
target = Target(frame, count).with_name("water")
```

Build a frame in memory (no PDB file) with `molrs.Frame.from_dict`:

```python
import molrs
import numpy as np

frame = molrs.Frame.from_dict({
    "blocks": {
        "atoms": {
            "x": np.array([0.00,  0.96, -0.24]),
            "y": np.array([0.00,  0.00,  0.93]),
            "z": np.zeros(3),
            "element": ["O", "H", "H"],
        }
    },
    "meta": {},
})
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

## Relaxation-assisted packing

Attach an in-loop relaxer to reshape a flexible molecule's reference
geometry *during* packing — useful for long chains that must fold to
fit. `with_relaxer` requires `count == 1` (every copy shares the
reference geometry the relaxer rewrites):

```python
from molpack import Target, TorsionMcRelaxer

chain = TorsionMcRelaxer(frame).with_steps(20).with_self_avoidance(1.5)
target = Target(frame, count=1).with_relaxer(chain)
```

Two relaxers ship built in: `TorsionMcRelaxer` (engine-free Monte-Carlo
torsion sampling, always available) and `LBFGSRelaxer` (force-field
L-BFGS, `ff` feature). See
[In-loop relaxers](../api-reference.md#in-loop-relaxers) for their
options.

## Per-target solver budget

Override the maximum perturbation budget for this target:

```python
target = target.with_perturb_budget(50)  # default: derived from count
```

Useful when one species is significantly harder to place than the rest.
