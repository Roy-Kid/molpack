# Restraints

Restraints are geometric regions (or half-spaces) that every atom of a
target — or a chosen subset — must lie inside. molpack ships five
built-in restraints.

## Built-ins

| Class                    | Constructor arguments                          | Meaning |
|--------------------------|------------------------------------------------|---------|
| `InsideBoxRestraint`     | `min: [x,y,z]`, `max: [x,y,z]`, `periodic=(False, False, False)` | axis-aligned box |
| `InsideSphereRestraint`  | `center: [x,y,z]`, `radius: float`             | closed ball |
| `OutsideSphereRestraint` | `center: [x,y,z]`, `radius: float`             | complement of closed ball |
| `AbovePlaneRestraint`    | `normal: [nx,ny,nz]`, `distance: float`        | half-space $\mathbf{n}\cdot\mathbf{x} \ge d$ |
| `BelowPlaneRestraint`    | `normal: [nx,ny,nz]`, `distance: float`        | half-space $\mathbf{n}\cdot\mathbf{x} \le d$ |

All arguments are standard Python floats / lists.

```python
from molpack import (
    AbovePlaneRestraint,
    BelowPlaneRestraint,
    InsideBoxRestraint,
    InsideSphereRestraint,
    OutsideSphereRestraint,
)

box   = InsideBoxRestraint([0, 0, 0], [40, 40, 40])
ball  = InsideSphereRestraint([0, 0, 0], 20.0)
shell = OutsideSphereRestraint([0, 0, 0], 10.0)
above = AbovePlaneRestraint(normal=[0, 0, 1], distance=5.0)
below = BelowPlaneRestraint(normal=[0, 0, 1], distance=20.0)
```

## Periodic boxes

`InsideBoxRestraint` doubles as the PBC declaration. Passing a
`periodic` tuple turns any subset of axes periodic:

```python
InsideBoxRestraint([0, 0, 0], [30, 30, 30], periodic=(True, True, True))
```

See [Periodic boundaries](periodic-boundaries.md) for the full
semantics and validation rules.

## Stacking multiple restraints

Apply several restraints to the same target by chaining `.with_restraint()`:

```python
target = (
    Target(frame, count=500)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
    .with_restraint(OutsideSphereRestraint([20, 20, 20], 5.0))
)
```

Each call attaches an independent restraint. All active restraints are
evaluated at every optimizer step.

## Scopes

A restraint can be applied at two scopes:

- **Whole target** — `target.with_restraint(r)` — penalises every atom.
- **Atom subset** — `target.with_atom_restraint([0, 1, 2], r)` —
  penalises only the listed atoms (0-based indices).

Example — lipid bilayer: pin heads above z=12, tails below z=2:

```python
lipid = (
    Target(frame, count=20)
    .with_name("lipid")
    .with_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 14]))
    .with_atom_restraint([0, 1],   AbovePlaneRestraint([0, 0, 1], 12.0))
    .with_atom_restraint([30, 31], BelowPlaneRestraint([0, 0, 1], 2.0))
)
```

## Global restraints

To apply one restraint to every target in a pack, attach it on the
packer:

```python
packer = (
    Molpack()
    .with_global_restraint(InsideBoxRestraint([0, 0, 0], [40, 40, 40]))
)
```

Semantically equivalent to calling `.with_restraint(r)` on every
target, but avoids the duplication.

## Custom restraints

Pass any object implementing `f(x, scale, scale2) -> float` and
`fg(x, scale, scale2) -> (float, (gx, gy, gz))`. See the
`Restraint` Protocol in `molpack` for the full contract.

```python
class SphereRestraint:
    def __init__(self, center, radius):
        self.c = np.asarray(center)
        self.r = radius
    def f(self, x, scale, scale2):
        d = np.linalg.norm(np.asarray(x) - self.c) - self.r
        return scale2 * d * d if d > 0 else 0.0
    def fg(self, x, scale, scale2):
        rel = np.asarray(x) - self.c
        d = float(np.linalg.norm(rel))
        over = d - self.r
        if over <= 0:
            return 0.0, (0.0, 0.0, 0.0)
        factor = 2 * scale2 * over / d
        return scale2 * over * over, tuple(factor * rel)

target = Target(frame, count=10).with_restraint(SphereRestraint([0,0,0], 6.0))
```

## Semantics

Every restraint contributes a continuously differentiable penalty
$f_{\text{rest}}(\mathbf{x})$ that is zero inside the allowed region
and rises quadratically outside. The aggregate objective minimised by
the packer is:

$$
U(\mathbf{x}) = f_{\text{dist}}(\mathbf{x}) + f_{\text{rest}}(\mathbf{x})
$$

where $f_{\text{dist}}$ is the pairwise distance-violation sum for the
user-specified `tolerance`. Convergence is declared when both fall
below `precision`.

!!! note "Restraints vs hard constraints"
    The five built-in classes are *soft penalties* — the optimizer may
    momentarily produce a violating configuration while searching. Hard
    geometric constraints (frozen placement, rotation bounds) are set
    on the `Target` directly via `fixed_at` and `with_rotation_bound`.
