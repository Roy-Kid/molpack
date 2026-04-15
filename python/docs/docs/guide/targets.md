# Targets

A `Target` describes one type of molecule to pack: its template
coordinates, van der Waals radii, element symbols, and the number of
copies to produce.

## Construction

```python
from molpack import Target

target = Target.from_coords(positions, radii, count, elements=None)
```

- `positions` — `(N, 3)` float64 numpy array of template coordinates.
- `radii`     — `(N,)` float64 numpy array, same length as positions.
- `count`     — number of copies to produce.
- `elements`  — optional list of element symbols, length $N$. Defaults
                to `"X"` per atom.

## Metadata

```python
target.with_name("water")  # labels the species in logs and traces
```

Read-only properties:

```python
target.natoms      # number of template atoms
target.count       # requested copies
target.elements    # list[str]
target.is_fixed    # True if placement is frozen (see below)
```

All builder methods are **immutable** — they return a new `Target` and
leave the original unchanged.

## Centering

By default the template is re-centered on its mean position before
packing (Packmol convention). Control explicitly:

```python
target = target.with_center()          # force centering on
target = target.without_centering()    # keep template origin
```

## Fixed placement

Pin a target at a specific location (useful for reference molecules
like a protein that should not move):

```python
target = target.fixed_at([10.0, 20.0, 30.0])                      # translate
target = target.fixed_at_with_euler(                              # + Euler
    position=[10.0, 20.0, 30.0],
    euler=[0.0, 1.57, 0.0],                                       # radians
)
```

Fixed targets contribute their atoms to distance-exclusion against
other species but are not themselves moved by the optimizer.

## Rotation bounds

For non-fixed targets you can restrict the rotational search window
about each axis (Packmol's `constrain_rotation` directive). Arguments
are in **degrees**:

```python
target = target.constrain_rotation_x(center_deg=0.0,  half_width_deg=15.0)
target = target.constrain_rotation_y(center_deg=90.0, half_width_deg=10.0)
target = target.constrain_rotation_z(center_deg=0.0,  half_width_deg=5.0)
```

## Attaching restraints

Two scopes:

### All atoms of the target

```python
from molpack import InsideBox

target = target.with_constraint(InsideBox([0, 0, 0], [40, 40, 40]))
```

### A subset of atoms

```python
from molpack import BelowPlane, AbovePlane

target = target.with_constraint_for_atoms(
    [31, 32],                                    # 1-based atom indices
    BelowPlane([0.0, 0.0, 1.0], 2.0),
)
```

!!! note "1-based indexing"
    `with_constraint_for_atoms` uses Packmol's 1-based atom indices.
    Passing `0` or an index `> natoms` raises `ValueError`.

Multiple restraints compose with `.and_()` (see [Restraints](restraints.md)).

## Per-target solver budget

Override the maximum optimizer iterations for this target:

```python
target = target.with_maxmove(50)  # default: derived from molecule size
```

Useful when one species is significantly harder to place than the rest.
