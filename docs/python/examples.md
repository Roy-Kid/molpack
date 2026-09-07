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
python examples/pack_peo_stl.py 25 200 130 42
```

Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the per-iteration progress log.
Open-space PEO defaults `LatticeGrow` then `GenCanPack.with_restart` at 2.0 Å.

`pack_peo_stl.py` runs `LatticeGrow` alone. Its cavity is the shipped
`examples/pack_peo/dendrite.stl` — a watertight dendrite, a trunk that forks
three times into 29 branches, 8 732 triangles, authored for a 130 Å cell with a
≈362 000 Å³ cavity — loaded with `scale` so one mesh serves any cell size. The
default 200 × EO25 (35 600 atoms) fills it at 1.026 g/cm³ in 2.2 s.

Three things that scene makes concrete:

- **The push-off is the wrong follow-up at melt density in a cavity.** It has
  nowhere to put the overlap it resolves except through the wall: on this scene
  `GenCanPack.with_restart` spent 1 h 45 min moving `fdist` 3.99 → 3.28 while
  `frest` went 0.42 → 6.43 — worst excursion 6.5 Å → 25 Å. Residual contacts at
  melt density are honest, and the force field downstream removes them.
- **The mask confines the walk, not the decorated atoms.** Decoration rebuilds
  the molecule from its own bonds and angles along the track and drifts off it;
  `with_track_tweak` is the lever. At EO25 the worst excursion after the walk is
  ≈36 Å at `0.0`, ≈16 Å at the `0.35` default and ≈3 Å at `1.5`, which is what
  the example asks for. 10% of atoms still end up to 6.5 Å outside, so a mesh
  that has to hold a wall is authored with the clearance already in it.
- **`precision` is a distance in disguise.** `frest` is the largest per-atom
  `0.01 · d²`, so `frest < precision` means `d < 10·√precision`: the default
  `1e-2` calls a run converged with an atom 1 Å outside a region.

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
