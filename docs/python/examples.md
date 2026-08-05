# Examples

Five canonical Packmol workloads ported to Python. Each lives under
`python/examples/` in the repo and is regression-tested against the
equivalent Rust example (same RNG seed → identical final coordinates).

| Script                | Packmol analogue  | What it shows |
|-----------------------|-------------------|---------------|
| `pack_water_cube.py`  | —                 | hello-world: 100 waters in a box, frame via `molrs.Frame.from_dict` |
| `pack_mixture.py`     | `mixture.inp`     | two species co-packed in one box |
| `pack_bilayer.py`     | `bilayer.inp`     | atom-subset restraints for layer-molecule orientation |
| `pack_interface.py`   | `interface.inp`   | fixed reference molecule + two solvents |
| `pack_spherical.py`   | `spherical.inp`   | nested spheres, double-layer shell |
| `pack_solvprotein.py` | `solvprotein.inp` | fixed solute solvated by water + ions |

Install molpack once; the `molrs` dependency comes with it:

```bash
pip install molcrafts-molpack
```

`pack_water_cube.py` builds its frame in memory with `molrs.Frame.from_dict`
(no PDB file); the others load PDB files via `molrs.io.read_pdb`.

## Running

```bash
cd molpack/python
pip install -e .
python examples/pack_water_cube.py       # no PDB file
python examples/pack_mixture.py          # requires molrs
```

Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the per-iteration progress log.

## Example: mixture

The `pack_mixture.py` example reproduces Packmol's classic `mixture.inp`:

```python
import molrs
from molpack import InsideBoxRestraint, Molpack, Target

water_frame = molrs.io.read_pdb("water.pdb")
urea_frame  = molrs.io.read_pdb("urea.pdb")

box = InsideBoxRestraint([0, 0, 0], [40, 40, 40])

water = Target(water_frame, count=1000).with_name("water").with_restraint(box)
urea  = Target(urea_frame,  count=400).with_name("urea").with_restraint(box)

packer = Molpack().with_tolerance(2.0).with_seed(1_234_567)
result = packer.pack_with_report([water, urea], max_loops=400)
print(f"converged={result.converged}  natoms={result.natoms}")
```

## Example: water cube

```python
import molrs
import numpy as np
from molpack import InsideBoxRestraint, Molpack, Target

frame = molrs.Frame.from_dict({
    "blocks": {
        "atoms": {
            "x": np.array([0.00,  0.9572, -0.2400]),
            "y": np.array([0.00,  0.0000,  0.9266]),
            "z": np.zeros(3),
            "element": ["O", "H", "H"],
        }
    },
    "meta": {},
})

water = Target(frame, count=100).with_name("water").with_restraint(
    InsideBoxRestraint([0, 0, 0], [30, 30, 30])
)
packer = Molpack().with_tolerance(2.0).with_log_level("quiet").with_seed(42)
result = packer.pack_with_report([water], max_loops=200)
print(f"converged={result.converged}  natoms={result.natoms}")
```
