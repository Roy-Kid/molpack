# Examples

Five canonical Packmol workloads ported to Python. Each lives under
`python/examples/` in the repo and is regression-tested against the
equivalent Rust example (same RNG seed → identical final coordinates).

| Script                  | Packmol analogue | What it shows |
|-------------------------|------------------|---------------|
| `pack_water_cube.py`    | —                | hello-world: 100 waters in a box |
| `pack_mixture.py`       | `mixture.inp`    | two species co-packed in one box |
| `pack_bilayer.py`       | `bilayer.inp`    | atom-subset restraints for lipid orientation |
| `pack_interface.py`     | `interface.inp`  | fixed reference molecule + two solvents |
| `pack_spherical.py`     | `spherical.inp`  | nested spheres, double-layer vesicle |
| `pack_solvprotein.py`   | `solvprotein.inp`| fixed protein solvated by water + ions |

All examples load PDB files via `molcrafts-molrs` (`molrs.read_pdb`)
and feed the coordinates + elements into `Target.from_coords`. If you
prefer not to install `molrs`, replace the loader with your own PDB
parser or plain numpy — `molpack` is I/O-free.

## Running

```bash
cd molpack/python
pip install -e .
pip install molcrafts-molrs          # for PDB I/O
python examples/pack_mixture.py
```

Set `MOLPACK_EXAMPLE_XYZ=1` to dump an XYZ trajectory under
`examples/out/`. Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the
per-iteration progress log.

## Example: mixture

The `pack_mixture.py` example reproduces Packmol's classic
`mixture.inp`:

```python
import numpy as np
import molrs
from molpack import Target, Packer, InsideBox

# Load PDB templates via molrs.
water = molrs.read_pdb("water.pdb")
urea  = molrs.read_pdb("urea.pdb")

# Extract (positions, elements) into numpy arrays.
def frame_to_arrays(frame):
    atoms = frame.get_block("atoms")
    pos = np.stack([atoms.get_float("x"),
                    atoms.get_float("y"),
                    atoms.get_float("z")], axis=1)
    els = atoms.get_string("element")
    rad = np.array([molrs.Box.vdw_radius(e) for e in els])  # or hardcode
    return pos, rad, els

water_pos, water_rad, water_els = frame_to_arrays(water)
urea_pos,  urea_rad,  urea_els  = frame_to_arrays(urea)

box = InsideBox([0, 0, 0], [40, 40, 40])

water_target = (
    Target.from_coords(water_pos, water_rad, count=1000, elements=water_els)
    .with_name("water")
    .with_constraint(box)
)
urea_target = (
    Target.from_coords(urea_pos, urea_rad, count=400, elements=urea_els)
    .with_name("urea")
    .with_constraint(box)
)

result = Packer().pack([water_target, urea_target], max_loops=400, seed=1_234_567)
print(f"converged={result.converged}  natoms={result.natoms}")
```

See the other example scripts in the repository for the full set.
