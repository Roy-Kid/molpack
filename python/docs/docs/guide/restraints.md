# Restraints

Restraints are geometric regions (or half-spaces) that every atom of a
target — or a chosen subset — must lie inside. molpack ships five
built-in restraints; any number can be composed.

## Built-ins

| Class           | Constructor arguments                         | Meaning |
|-----------------|-----------------------------------------------|---------|
| `InsideBox`     | `min: [x, y, z]`, `max: [x, y, z]`            | axis-aligned box |
| `InsideSphere`  | `radius: float`, `center: [x, y, z]`          | closed ball |
| `OutsideSphere` | `radius: float`, `center: [x, y, z]`          | complement of closed ball |
| `AbovePlane`    | `normal: [nx, ny, nz]`, `distance: float`     | half-space $\mathbf{n}\cdot\mathbf{x} \ge d$ |
| `BelowPlane`    | `normal: [nx, ny, nz]`, `distance: float`     | half-space $\mathbf{n}\cdot\mathbf{x} \le d$ |

All arguments are standard Python floats / lists. Internally they are
forwarded to the corresponding Rust `…Restraint` struct.

```python
from molpack import InsideBox, InsideSphere, OutsideSphere, AbovePlane, BelowPlane

box      = InsideBox([0, 0, 0], [40, 40, 40])
ball     = InsideSphere(radius=20.0, center=[0, 0, 0])
shell    = OutsideSphere(radius=10.0, center=[0, 0, 0])
above    = AbovePlane(normal=[0, 0, 1], distance=5.0)
below    = BelowPlane(normal=[0, 0, 1], distance=20.0)
```

## Composition

Combine restraints with `.and_()` to form a `MoleculeConstraint` — a
bundle that is attached atomically:

```python
from molpack import InsideBox, OutsideSphere

bundle = InsideBox([0, 0, 0], [40, 40, 40]).and_(
    OutsideSphere(radius=5.0, center=[20, 20, 20])
)

target = target.with_constraint(bundle)
```

Chain further:

```python
from molpack import InsideSphere, AbovePlane, BelowPlane

shell_band = (
    InsideSphere(40.0, [0, 0, 0])
    .and_(AbovePlane([0, 0, 1], 10.0))
    .and_(BelowPlane([0, 0, 1], 30.0))
)
```

Calling `and_` on a `MoleculeConstraint` (the composite) also works
and produces a new composite.

## Scopes

A restraint can be applied at three scopes. See [Targets](targets.md)
for the detailed per-target API and the `Packer` guide for broadcast.

- **Whole target** — `target.with_constraint(r)`.
- **Atom subset** — `target.with_constraint_for_atoms([1, 2, 3], r)`
  (1-based indices, Packmol convention).
- **All targets** — `packer.pack(…)` takes no broadcast argument in
  the Python binding; replicate by calling `with_constraint` on each
  target.

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

!!! tip "Why _restraints_ and not _constraints_?"
    In the Rust crate each of the five built-ins is a soft penalty,
    not a hard constraint — the optimizer may momentarily produce a
    violating configuration while searching. The umbrella Python class
    `MoleculeConstraint` is named after Packmol's historical
    "constraint" keyword but its semantics are those of a restraint
    bundle.
