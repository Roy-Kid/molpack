# Restraints

Restraints are geometric regions (or half-spaces) that every atom of a
target — or a chosen subset — must lie inside. molpack ships five
built-in restraints.

## Built-ins

| Class           | Constructor arguments                     | Meaning |
|-----------------|-------------------------------------------|---------|
| `InsideBox`     | `min: [x,y,z]`, `max: [x,y,z]`           | axis-aligned box |
| `InsideSphere`  | `radius: float`, `center: [x,y,z]`       | closed ball |
| `OutsideSphere` | `radius: float`, `center: [x,y,z]`       | complement of closed ball |
| `AbovePlane`    | `normal: [nx,ny,nz]`, `distance: float`  | half-space $\mathbf{n}\cdot\mathbf{x} \ge d$ |
| `BelowPlane`    | `normal: [nx,ny,nz]`, `distance: float`  | half-space $\mathbf{n}\cdot\mathbf{x} \le d$ |

All arguments are standard Python floats / lists.

```python
from molpack import InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane

box   = InsideBox([0, 0, 0], [40, 40, 40])
ball  = InsideSphere(radius=20.0, center=[0, 0, 0])
shell = OutsideSphere(radius=10.0, center=[0, 0, 0])
above = AbovePlane(normal=[0, 0, 1], distance=5.0)
below = BelowPlane(normal=[0, 0, 1], distance=20.0)
```

## Stacking multiple restraints

Apply several restraints to the same target by chaining `.with_restraint()`:

```python
target = (
    Target("water", frame, count=500)
    .with_restraint(InsideBox([0, 0, 0], [40, 40, 40]))
    .with_restraint(OutsideSphere(5.0, [20, 20, 20]))
)
```

Each call attaches an independent restraint. All active restraints are
evaluated at every optimizer step.

## Scopes

A restraint can be applied at two scopes:

- **Whole target** — `target.with_restraint(r)` — penalises every atom.
- **Atom subset** — `target.with_restraint_for_atoms([1, 2, 3], r)` —
  penalises only the listed atoms (1-based Packmol indices).

Example — lipid bilayer: pin heads above z=12, tails below z=2:

```python
lipid = (
    Target("lipid", frame, count=20)
    .with_restraint(InsideBox([0, 0, 0], [40, 40, 14]))
    .with_restraint_for_atoms([1, 2],   AbovePlane([0, 0, 1], 12.0))
    .with_restraint_for_atoms([31, 32], BelowPlane([0, 0, 1], 2.0))
)
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

!!! note "Restraints vs constraints"
    The five built-in classes are *soft penalties* — the optimizer may
    momentarily produce a violating configuration while searching. Hard
    geometric constraints (frozen placement, rotation bounds) are set
    on the `Target` directly via `fixed_at` and `constrain_rotation_*`.
