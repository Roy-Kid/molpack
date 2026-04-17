# Getting Started

A minimal end-to-end pack: 100 water molecules inside a 40 Å cube.

## 1. Load a molecule

Use `molrs.read_pdb` to load a template PDB file — the returned
`Frame` can be passed directly to `Target`:

```python
import molrs

frame = molrs.read_pdb("water.pdb")
```

No `molrs`? Build the frame as a plain dict:

```python
import numpy as np

frame = {
    "atoms": {
        "x": np.array([0.00,  0.96, -0.24]),
        "y": np.array([0.00,  0.00,  0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
}
```

## 2. Create a Target

A `Target` bundles a molecule template with the number of copies to pack.
VdW radii are looked up automatically from element symbols (Bondi 1964).

```python
from molpack import Target

water = Target(frame, count=100).with_name("water")
```

Arguments:

- `frame` — any object supporting `frame["atoms"]`, with columns
  `"x"`, `"y"`, `"z"`, and `"element"` (or `"symbol"` for molrs PDB frames).
- `count` — number of copies to produce.

A display label is optional — attach one via `.with_name("...")`.

All builder methods are **immutable** — they return a new `Target`.

## 3. Attach a restraint

Every target needs at least one restraint — the geometric region it
should be packed into.

```python
from molpack import InsideBoxRestraint

water = water.with_restraint(
    InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])
)
```

Five built-in restraints: `InsideBoxRestraint`, `InsideSphereRestraint`,
`OutsideSphereRestraint`, `AbovePlaneRestraint`, `BelowPlaneRestraint`.
Stack multiple restraints with repeated `.with_restraint()` calls —
see [Restraints](guide/restraints.md).

## 4. Pack

```python
from molpack import Molpack

packer = Molpack().with_tolerance(2.0).with_seed(42)
result = packer.pack([water], max_loops=200)

print(f"converged: {result.converged}")
print(f"natoms:    {result.natoms}")
print(f"fdist:     {result.fdist:.4f}   (distance-violation sum)")
print(f"frest:     {result.frest:.4f}   (restraint-violation sum)")
```

`result.positions` is an `(N, 3)` float64 numpy array of packed
coordinates. `result.elements` is the element list in the same order.

## 5. Save

`molpack` does not write files directly — Frame is the canonical
output. Hand `result.frame` (or positions/elements) to a writer:

```python
import molrs

# result.frame is a dict compatible with molrs.Frame. Compose a
# writable frame and persist:
out = molrs.Frame()
atoms = molrs.Block()
f = result.frame["atoms"]
atoms.insert("x", f["x"])
atoms.insert("y", f["y"])
atoms.insert("z", f["z"])
atoms.insert("element", f["element"])
out["atoms"] = atoms
molrs.write_xyz("packed.xyz", out)
```

Or use numpy directly (`np.savetxt`, `.npz`, …) on `result.positions`.

## Full script

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

result = (
    Molpack().with_tolerance(2.0).with_seed(42).pack([water], max_loops=200)
)

assert result.converged
print(f"packed {result.natoms} atoms")
```

## Next steps

- [Targets](guide/targets.md) — orientation, centering, fixed placement.
- [Restraints](guide/restraints.md) — per-atom scoping, stacking.
- [Packer](guide/packer.md) — all builder options.
- [Examples](examples.md) — five complete Packmol workloads.
