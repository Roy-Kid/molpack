# Python Periodic Boundaries

Periodic boundary conditions in the Python binding come in two flavours
and you can use whichever one reads more naturally for your setup:

- **Packer-level** — `Molpack.with_periodic_box(min, max)`. Every axis
  is periodic. This is the Packmol `pbc` keyword; it is also what
  `load_script` wires up for you automatically when an `.inp` script
  contains a `pbc` directive.
- **Restraint-level** — `InsideBoxRestraint(min, max, periodic=(...))`.
  Per-axis control, so you can mark only some axes as wrapping (e.g.
  slab geometries).

Only one periodic box may exist per run. If both paths are used they
must agree exactly (same bounds, same per-axis flags) or `pack()`
raises `ConflictingPeriodicBoxesError`.

## A Fully Periodic Box, on the packer

```python
from molpack import Molpack

packer = Molpack().with_periodic_box(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 30.0],
)
```

That is the right choice when the box is a global property of the run
(all molecules share the same domain) — which is the common case, and
the only one `.inp` scripts can express.

## A Fully Periodic Box, on an `InsideBoxRestraint`

```python
from molpack import InsideBoxRestraint

box = InsideBoxRestraint(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 30.0],
    periodic=(True, True, True),
)
```

That says two things at once:

1. atoms are softly confined to the box;
2. pair distances wrap through that box on the periodic axes.

## Slabs and Mixed Boundary Conditions

You can also switch periodicity on axis by axis:

```python
slab = InsideBoxRestraint(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 100.0],
    periodic=(True, True, False),
)
```

That is the usual setup for slab geometries: periodic in-plane, open along the
surface normal.

## How `molpack` Interprets PBC

At `pack()` time, the solver gathers periodic-box declarations from
two sources:

1. the packer itself (`with_periodic_box` or a script `pbc` directive);
2. every restraint on every target (`InsideBoxRestraint` with any axis
   marked `periodic`).

Then:

- If there are none, the run is non-periodic.
- If there is exactly one box declaration (from either source), that
  box defines the system cell.
- If several declarations exist, they must agree exactly on bounds
  and per-axis flags.

Mismatched periodic boxes raise `ConflictingPeriodicBoxesError`.

## From a `.inp` Script

Scripts use Packmol's `pbc` keyword, with two accepted forms:

```text
pbc 30.0 30.0 30.0                    # min = [0, 0, 0], max = [30, 30, 30]
pbc -15.0 -15.0 -15.0  15.0 15.0 15.0 # explicit min / max
```

`load_script` turns this into a packer-level periodic box — the same
thing `Molpack.with_periodic_box(min, max)` would have produced by
hand. Without a `pbc` or an `inside` restraint the packer has no
spatial constraint and falls back to inferring a box from the initial
random placement, which drives the cell grid to tens of millions of
cells; always give the packer a spatial constraint.

## A Box Must Be Valid

Boxes with non-positive extents are rejected with `InvalidPBCBoxError`. This is
usually a sign that the min and max corners were swapped or one axis was
accidentally collapsed.

## A Typical Pattern

For most MD-style setups, the pattern is simply:

```python
from molpack import Molpack, Target

target = target.with_restraint(
    InsideBoxRestraint(
        [0.0, 0.0, 0.0],
        [30.0, 30.0, 30.0],
        periodic=(True, True, True),
    )
)

result = Molpack().with_seed(42).pack([target], max_loops=200)
```

## Takeaway

Pick the path that matches the shape of your problem:

- one global box shared by everyone → `Molpack.with_periodic_box` (or
  a `pbc` line in the `.inp` script);
- per-axis control, slab geometries, or per-target bounds → set
  `periodic=` on `InsideBoxRestraint`.

Either way, the packer will derive a single consistent system cell.
