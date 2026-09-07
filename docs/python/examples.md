# Examples

Packmol workloads ported to Python live under `python/examples/`. The five
`.inp` analogues are regression-tested against the equivalent Rust example
(same RNG seed → identical final coordinates). Polymer scenes sit beside
them: chemistry from molrs SMILES + molpy `PolymerBuilder`, packing from
molpack — no hand-placed coordinates.

| Script                | Packmol analogue  | What it shows |
|-----------------------|-------------------|---------------|
| `pack_water_cube.py`  | —                 | hello-world: 100 waters in a box, frame via `molrs.Frame` |
| `pack_mixture.py`     | `mixture.inp`     | two species co-packed in one box |
| `pack_bilayer.py`     | `bilayer.inp`     | atom-subset restraints for layer-molecule orientation |
| `pack_interface.py`   | `interface.inp`   | fixed reference molecule + two solvents |
| `pack_spherical.py`   | `spherical.inp`   | nested spheres, double-layer shell |
| `pack_solvprotein.py` | `solvprotein.inp` | fixed solute solvated by water + ions |
| `pack_peo_linear.py`  | —                 | open-space linear PEO: `LatticeGrow` @ 2.0 Å then `GenCanPack.with_restart` |
| `pack_peo_mix.py`     | —                 | linear + 4-arm star, two `Target`s, one box, one `LatticeGrow.run` |
| `pack_peo_topo.py`    | —                 | 4-arm star (`LatticeGrow`) and ring (named reject, then `GenCanPack`) |
| `pack_peo_stl.py`     | —                 | linear PEO inside a branched STL cavity (`StlRegion` masks `LatticeGrow` sites) |

Install molpack once; the `molrs` dependency comes with it:

```bash
pip install molcrafts-molpack
```

Each script is standalone: no shared helper. `pack_water_cube.py` builds
its frame in memory with `molrs.Frame` (no PDB file). The Packmol-port
scripts load PDB files via `molrs.io.read_pdb`. The `pack_peo_*.py`
scenes build polymers from SMILES + `PolymerBuilder` instead. Writes
go through molrs (`molrs.io.mrec.write_frame`, `write_lammps_traj`,
`write_lammps_dump_local`).

## Running

```bash
cd molpack/python
pip install -e .
python examples/pack_water_cube.py       # no PDB file
python examples/pack_mixture.py          # requires molrs
python examples/pack_peo_linear.py 8 8 0.5 42
python examples/pack_peo_mix.py 4 2 4 4 0.5 42
python examples/pack_peo_topo.py star 4 8 0.5 42
python examples/pack_peo_stl.py 4 4 48 42
```

Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the per-iteration progress log.
Open-space PEO defaults `LatticeGrow` then `GenCanPack.with_restart` at
2.0 Å; `pack_peo_stl.py` uses the same pipeline with `StlRegion` masking
diamond sites outside the mesh. Its cavity is the shipped
`examples/pack_peo/dendrite.stl` — a watertight dendrite, a trunk that
forks three times into 29 branches — loaded with `scale` so one mesh
serves any cell size. It pushes off at `precision=1e-4`: the default
1e-2 leaves an atom ~1 Å outside a region and still reports `converged`,
because `frest` is the largest per-atom `0.01 · d²`. The lattice mask
confines the *walk*, not the decorated atoms — a mesh that has to hold a
wall should be authored with the clearance already in it. The drift grows
with the backbone (≈1 Å at EO3, ≈4 Å at EO4, ≈8 Å at EO5); the push-off
still recovers EO4, and from EO5 up it reports `converged=False`.

Each example writes its outputs to `python/examples/out/` (created on
demand, git-ignored) — the path is script-relative, so the working
directory does not matter:

- `{stem}.mrec` — molrs scientific record (`molrs.io.mrec.write_frame`)
- `{stem}.lammpstrj` — LAMMPS dump custom (OVITO particle topology)
- `{stem}.dump.local` — LAMMPS dump local bonds (`batom1`/`batom2`), for
  OVITO [Load trajectory](https://www.ovito.org/manual/reference/pipelines/modifiers/load_trajectory.html)
  (open the `.lammpstrj`, then overlay the dump local file). Skipped when
  the packed frame has no bonds (the water-cube template).

`StlRegion` answers the two region questions for a batch of points, so a
caller can check what the packer was told to enforce:

```python
cavity = molpack.StlRegion.from_file("dendrite.stl")
inside = cavity.contains(state.positions)       # (n,) bool
depth = cavity.signed_distance(state.positions)  # (n,) Å, negative inside
```

## Example: mixture

The `pack_mixture.py` example reproduces Packmol's classic `mixture.inp`:

```python
import molrs
from molpack import GenCanPack, InsideBoxRestraint, Target

water_frame = molrs.io.read_pdb("water.pdb")
urea_frame  = molrs.io.read_pdb("urea.pdb")

box = InsideBoxRestraint([0, 0, 0], [40, 40, 40])

water = Target(water_frame, count=1000).with_name("water").with_restraint(box)
urea  = Target(urea_frame,  count=400).with_name("urea").with_restraint(box)

packer = GenCanPack().with_tolerance(2.0).with_seed(1_234_567)
result = packer.run([water, urea], max_loops=400)
print(f"converged={result.converged}  natoms={result.natoms}")
```

## Example: water cube

```python
import molrs
import numpy as np
from molpack import GenCanPack, InsideBoxRestraint, Target

frame = molrs.Frame({
    "atoms": {
        "x": np.array([0.00,  0.9572, -0.2400]),
        "y": np.array([0.00,  0.0000,  0.9266]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
})

water = Target(frame, count=100).with_name("water").with_restraint(
    InsideBoxRestraint([0, 0, 0], [30, 30, 30])
)
packer = GenCanPack().with_tolerance(2.0).with_progress(False).with_seed(42)
result = packer.run([water], max_loops=200)
print(f"converged={result.converged}  natoms={result.natoms}")
```
