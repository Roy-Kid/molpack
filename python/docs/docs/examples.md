# Examples

Five canonical Packmol workloads ported to Python. Each lives under
`python/examples/` in the repo and is regression-tested against the
equivalent Rust example (same RNG seed → identical final coordinates).

| Script                | Packmol analogue  | What it shows |
|-----------------------|-------------------|---------------|
| `pack_water_cube.py`  | —                 | hello-world: 100 waters in a box, no I/O library |
| `pack_mixture.py`     | `mixture.inp`     | two species co-packed in one box |
| `pack_bilayer.py`     | `bilayer.inp`     | atom-subset restraints for lipid orientation |
| `pack_interface.py`   | `interface.inp`   | fixed reference molecule + two solvents |
| `pack_spherical.py`   | `spherical.inp`   | nested spheres, double-layer vesicle |
| `pack_solvprotein.py` | `solvprotein.inp` | fixed protein solvated by water + ions |

Examples that load PDB files (`pack_mixture.py` and above) require
`molcrafts-molrs`:

```bash
pip install molcrafts-molrs
```

`pack_water_cube.py` is fully standalone — it constructs the frame as a
plain dict.

## Running

```bash
cd molpack/python
pip install -e .
python examples/pack_water_cube.py       # standalone
python examples/pack_mixture.py          # requires molrs
```

Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the per-iteration progress log.

## Example: mixture

The `pack_mixture.py` example reproduces Packmol's classic `mixture.inp`:

```python
import molrs
from molpack import InsideBox, Molpack, Target

water_frame = molrs.read_pdb("water.pdb")
urea_frame  = molrs.read_pdb("urea.pdb")

box = InsideBox([0, 0, 0], [40, 40, 40])

water = Target("water", water_frame, count=1000).with_restraint(box)
urea  = Target("urea",  urea_frame,  count=400).with_restraint(box)

result = Molpack().pack([water, urea], max_loops=400, seed=1_234_567)
print(f"converged={result.converged}  natoms={result.natoms}")
```

## Example: water cube (standalone)

```python
import numpy as np
from molpack import InsideBox, Molpack, Target

frame = {
    "atoms": {
        "x": np.array([0.00,  0.9572, -0.2400]),
        "y": np.array([0.00,  0.0000,  0.9266]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
}

water = Target("water", frame, count=100).with_restraint(
    InsideBox([0, 0, 0], [30, 30, 30])
)
result = Molpack(tolerance=2.0, progress=False).pack([water], max_loops=200, seed=42)
print(f"converged={result.converged}  natoms={result.natoms}")
```
