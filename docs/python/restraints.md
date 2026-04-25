# Python Restraints

Restraints answer the question "where is this molecule allowed to go?" In
practice, they are the part of the API that makes a packing job feel like your
actual physical setup rather than a generic optimization problem.

## The Built-In Restraints

`molpack` ships five built-in restraint classes in Python:

| Class | Meaning |
|---|---|
| `InsideBoxRestraint(min, max, periodic=...)` | axis-aligned box |
| `InsideSphereRestraint(center, radius)` | closed ball |
| `OutsideSphereRestraint(center, radius)` | everything outside a sphere |
| `AbovePlaneRestraint(normal, distance)` | half-space above a plane |
| `BelowPlaneRestraint(normal, distance)` | half-space below a plane |

```python
from molpack import (
    AbovePlaneRestraint,
    BelowPlaneRestraint,
    InsideBoxRestraint,
    InsideSphereRestraint,
    OutsideSphereRestraint,
)
```

## Whole-Target vs Atom-Subset Scope

There are two common ways to attach a restraint:

```python
target = target.with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
target = target.with_atom_restraint([0, 1], AbovePlaneRestraint([0, 0, 1], 12.0))
```

Use whole-target scope when every atom should obey the same region. Use
atom-subset scope when only selected atoms matter, such as lipid head groups in
bilayer setups.

One easy gotcha: atom subsets use `0`-based indices, even if you are porting a
Packmol input deck that used `1`-based numbering.

## Stacking Restraints

Real systems usually need more than one geometric rule. Chaining is fine:

```python
target = (
    target
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
    .with_restraint(OutsideSphereRestraint([20, 20, 20], 5.0))
)
```

The solver treats restraints as soft penalties. In other words, stacking them
adds geometric pressure; it does not create a magical hard intersection outside
the objective function.

## Global Restraints

If every target should share the same boundary, attach it to the packer once:

```python
from molpack import InsideBoxRestraint, Molpack

packer = Molpack().with_global_restraint(
    InsideBoxRestraint([0, 0, 0], [40, 40, 40])
)
```

This is equivalent to attaching the same restraint to each target, just less
repetitive.

## Custom Restraints

You can also pass any Python object that provides:

- `f(x, scale, scale2) -> float`
- `fg(x, scale, scale2) -> (float, (gx, gy, gz))`

That makes quick experiments easy:

```python
import numpy as np

class SphereCap:
    def __init__(self, center, radius):
        self.center = np.asarray(center, dtype=float)
        self.radius = float(radius)

    def f(self, x, scale, scale2):
        rel = np.asarray(x, dtype=float) - self.center
        over = np.linalg.norm(rel) - self.radius
        return scale2 * over * over if over > 0 else 0.0

    def fg(self, x, scale, scale2):
        rel = np.asarray(x, dtype=float) - self.center
        dist = float(np.linalg.norm(rel))
        over = dist - self.radius
        if over <= 0:
            return 0.0, (0.0, 0.0, 0.0)
        grad = 2.0 * scale2 * over * rel / dist
        return scale2 * over * over, tuple(grad)
```

If you go custom, validate the gradient early. A restraint with the wrong `fg`
method can look fine in code and still make the optimizer miserable.

## Periodic Boxes Live on `InsideBoxRestraint`

Periodic boundaries are declared on the box itself:

```python
box = InsideBoxRestraint(
    [0.0, 0.0, 0.0],
    [30.0, 30.0, 100.0],
    periodic=(True, True, False),
)
```

That rule is explained in more detail in [Periodic Boundaries](periodic_boundaries.md).

## Takeaway

Restraints are where the chemistry starts to feel real. Start with the built-in
shapes, use atom subsets when orientation matters, and only reach for custom
Python restraint objects when the built-ins stop matching the geometry you need.
