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
frame = packer.pack([water], max_loops=200)

print(frame["atoms"].nrows)
```

`pack()` returns a ready-to-use `molrs.Frame`. If you need structured
diagnostics, call `pack_with_report()` instead; it returns a
`PackResult` with `.converged`, `.fdist`, `.frest`, `.positions`, and
`.frame`.

## 5. Save

`molpack` does not write files directly — Frame is the canonical
output. Hand the returned frame to a writer:

```python
import molrs

molrs.write_xyz("packed.xyz", frame)
```

Or use `pack_with_report()` and write `result.frame` if you also need
the diagnostic fields.

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

frame = (
    Molpack().with_tolerance(2.0).with_seed(42).pack([water], max_loops=200)
)

print(f"packed {frame['atoms'].nrows} atoms")
```

## Next steps

- [Targets](guide/targets.md) — orientation, centering, fixed placement.
- [Restraints](guide/restraints.md) — per-atom scoping, stacking.
- [Packer](guide/packer.md) — all builder options.
- [Examples](examples.md) — five complete Packmol workloads.
