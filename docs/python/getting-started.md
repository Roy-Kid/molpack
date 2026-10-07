# Quickstart

A minimal end-to-end pack: 100 water molecules inside a 40 Å cube.

## 1. Load a molecule

Use `molrs.io.read_pdb` to load a template PDB file — the returned
`Frame` can be passed directly to `Target`:

```python
import molrs

frame = molrs.io.read_pdb("water.pdb")
```

No PDB file? Build a `molrs.core.Frame` from arrays:

```python
import molrs
import numpy as np

frame = molrs.core.Frame({
    "atoms": {
        "x": np.array([0.00,  0.96, -0.24]),
        "y": np.array([0.00,  0.00,  0.93]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
})
```

## 2. Create a Target

A `Target` bundles a molecule template with the number of copies to pack.
VdW radii are looked up automatically from element symbols (Bondi 1964).

```python
from molpack import Target

water = Target(frame, count=100).with_name("water")
```

Arguments:

- `frame` — a `molrs.core.Frame` (`molpy.Frame` is the same class) with columns
  `"x"`, `"y"`, `"z"`, and `"element"`.
- `count` — number of copies to produce.

A display label is optional — attach one via `.with_name("...")`.

All builder methods are **immutable** — they return a new `Target`.

## 3. Attach a restraint

Every target needs at least one restraint — the geometric region it
should be packed into.

```python
import molrs

water = water.with_restraint(
    molrs.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])   # origin, lengths
)
```

Any molrs region is a restraint: `Sphere`, `Cuboid`, `Parallelepiped`,
`HalfSpace`, `Cylinder`, `Ellipsoid`, `Polyhedron`, `SphereUnion`, and
their `&` / `|` / `~` compositions — plus a family of collective
distribution-matching restraints. Stack multiple restraints with
repeated `.with_restraint()` calls — see
[Restraints](guide/restraints.md).

## 4. Pack

```python
from molpack import GencanPack

packer = GencanPack().with_tolerance(2.0).with_seed(42)
result = packer.run([water], max_loops=200)
frame = result.frame

print(frame["atoms"].n_rows)
```

`GencanPack` is the rigid-body entry — you choose the packing algorithm
by choosing the entry, and `CbmcGrow` is the chain-growth one. Both have
the same builders and the same terminal verb, `run()`, which returns a
`State` with `.frame`, `.converged`, `.fdist`, `.frest`,
`.positions`, `.degraded`, and `.intra`. An entry runs once: build a
new one for each pack.

## 5. Save

`molpack` does not write files directly — Frame is the canonical
output. Hand the returned frame to a writer:

```python
import molrs

molrs.io.write_xyz("packed.xyz", frame)
```

`result.frame` is the same object — keep the `State` around when you
need the diagnostic fields, or to continue with
`GencanPack().with_restart(result)` / `Target.fixed_from(result)`.

## Full script

```python
import molrs
from molpack import GencanPack, Target

frame = molrs.io.read_pdb("water.pdb")

water = (
    Target(frame, count=100)
    .with_name("water")
    .with_restraint(molrs.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0]))
)

result = (
    GencanPack().with_tolerance(2.0).with_seed(42).run([water], max_loops=200)
)

print(f"packed {result.frame['atoms'].n_rows} atoms")
```

## Next steps

- [Targets](guide/targets.md) — orientation, centering, fixed placement.
- [Restraints](guide/restraints.md) — per-atom scoping, stacking.
- [Packer](guide/packer.md) — all builder options.
- [Examples](examples.md) — five complete Packmol workloads.
