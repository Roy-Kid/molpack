# Getting Started

A minimal end-to-end pack: 100 water molecules inside a 40 Å cube.

## 1. Build a Target

A `Target` bundles the template geometry of one molecule type with the
number of copies to pack.

```python
import numpy as np
from molpack import Target

water_positions = np.array([
    [0.00,  0.00, 0.00],   # O
    [0.96,  0.00, 0.00],   # H
    [-0.24, 0.93, 0.00],   # H
])
water_radii = np.array([1.52, 1.20, 1.20])

water = Target.from_coords(
    water_positions, water_radii,
    count=100,
    elements=["O", "H", "H"],
).with_name("water")
```

Arguments:

- `positions`  — `(N, 3)` float64 numpy array, Å.
- `radii`      — `(N,)` float64 numpy array, Å (typically van der Waals).
- `count`      — number of copies to produce.
- `elements`   — optional element symbols; default `"X"` per atom.

## 2. Attach a restraint

Every target needs at least one restraint — the geometric region into
which it should be packed.

```python
from molpack import InsideBox

water = water.with_constraint(
    InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])
)
```

Five concrete restraints are shipped: `InsideBox`, `InsideSphere`,
`OutsideSphere`, `AbovePlane`, `BelowPlane`. Compose multiple restraints
with `.and_()` — see [Restraints](guide/restraints.md).

## 3. Pack

```python
from molpack import Packer

packer = Packer(tolerance=2.0, precision=0.01)
result = packer.pack([water], max_loops=200, seed=42)

print(f"converged: {result.converged}")
print(f"natoms:    {result.natoms}")
print(f"fdist:     {result.fdist:.4f}   (distance-violation sum)")
print(f"frest:     {result.frest:.4f}   (restraint-violation sum)")
```

Arguments:

- `targets`    — list of `Target`s (one per molecule type).
- `max_loops`  — outer-iteration budget per phase.
- `seed`       — optional reproducibility seed.

`result.positions` is an `(N, 3)` float64 numpy array of the packed
atom coordinates. `result.elements` is the element list in the same
order.

## 4. Save

`molpack` does not write files directly. Either:

- Use numpy (`np.savetxt`, `.npz`, …) on `result.positions`, or
- Pass coordinates to `molcrafts-molrs` to write PDB / XYZ:

```python
import molrs

frame = molrs.Frame()
atoms = molrs.Block()
atoms.set_float("x", result.positions[:, 0])
atoms.set_float("y", result.positions[:, 1])
atoms.set_float("z", result.positions[:, 2])
atoms.set_string("element", result.elements)
frame.set_block("atoms", atoms)
molrs.write_xyz("packed.xyz", frame)
```

## Full script

```python
import numpy as np
from molpack import Target, Packer, InsideBox

water_positions = np.array([
    [0.00,  0.00, 0.00],
    [0.96,  0.00, 0.00],
    [-0.24, 0.93, 0.00],
])
water_radii = np.array([1.52, 1.20, 1.20])

water = (
    Target.from_coords(water_positions, water_radii, count=100,
                       elements=["O", "H", "H"])
    .with_name("water")
    .with_constraint(InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

result = Packer(tolerance=2.0, precision=0.01).pack(
    [water], max_loops=200, seed=42,
)

assert result.converged
print(f"packed {result.natoms} atoms")
```

## Next steps

- [Targets](guide/targets.md) — orientation, centering, fixed placement.
- [Restraints](guide/restraints.md) — composition, per-atom scoping.
- [Packer](guide/packer.md) — all builder options.
- [Examples](examples.md) — five complete Packmol workloads ported to Python.
