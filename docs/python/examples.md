# Examples

Packmol workloads ported to Python live under `python/examples/`. The five
`.inp` analogues are regression-tested against the equivalent Rust example
(same RNG seed → identical final coordinates). Polymer scenes sit beside
them: chemistry and architecture from molrs (CGsmiles units + conformer,
grown by `molrs.builder.Assembler`), packing from molpack — no hand-placed
coordinates.

| Script                | Packmol analogue  | What it shows |
|-----------------------|-------------------|---------------|
| `pack_water_cube.py`  | —                 | hello-world: 100 waters in a box, frame via `molrs.core.Frame` |
| `pack_mixture.py`     | `mixture.inp`     | two species co-packed in one box |
| `pack_bilayer.py`     | `bilayer.inp`     | atom-subset restraints for layer-molecule orientation |
| `pack_interface.py`   | `interface.inp`   | fixed reference molecule + two solvents |
| `pack_spherical.py`   | `spherical.inp`   | nested spheres, double-layer shell |
| `pack_solvprotein.py` | `solvprotein.inp` | fixed solute solvated by water + ions |
| `pack_ion_dispersion.py` | —              | `SelfSeparation`: stop one species clustering, with a no-restraint control |
| `pack_peo_linear.py`  | —                 | open-space linear PEO: `LatticeGrow` @ 2.0 Å then `GencanPack.with_restart` |
| `pack_peo_mix.py`     | —                 | linear + 4-arm star, two `Target`s, one box, one `LatticeGrow.run` |
| `pack_peo_topo.py`    | —                 | 4-arm star (`LatticeGrow`) and ring (named reject, then `GencanPack`) |
| `pack_peo_mesh.py`    | —                 | linear PEO inside a branched mesh cavity (a `molrs.core.Polyhedron` masks `LatticeGrow` sites) |
| `pack_peo_void.py`    | —                 | linear PEO through the solvent-accessible void of a bead-spring frame (`~molrs.core.SphereUnion`) |

Install molpack once; the `molrs` dependency comes with it, and every
script needs nothing else:

```bash
pip install molcrafts-molpack
```

Each script is standalone: no shared helper. `pack_water_cube.py` builds
its frame in memory with `molrs.core.Frame` (no PDB file). The Packmol-port
scripts load PDB files via `molrs.io.read_pdb`. The `pack_peo_*.py`
scenes build polymers from CGsmiles units
(`molrs.io.smiles.SmilesIr.from_fragment(body).to_template()`, given 3D
coordinates by `molrs.conformer.Conformer`) grown by
`molrs.builder.Assembler` with `molrs.builder.GrowthPlacer` instead. Writes
go through molrs (`molrs.io.write_mrec_frame`, `write_lammps_trajectory`,
`write_lammps_dump_local`).

## Running

```bash
cd molpack/python
pip install -e .
python examples/pack_water_cube.py       # no PDB file
python examples/pack_mixture.py          # requires molrs
python examples/pack_ion_dispersion.py   # packs twice: without, then with, the bound
python examples/pack_peo_linear.py 8 8 0.5 42
python examples/pack_peo_mix.py 4 2 4 4 0.5 42
python examples/pack_peo_topo.py star 4 8 0.5 42
python examples/pack_peo_mesh.py 25 200 130 42
python examples/pack_peo_void.py frame.data 25 200 42
```

Set `MOLPACK_EXAMPLE_PROGRESS=0` to suppress the per-iteration progress log.
Open-space PEO defaults `LatticeGrow` then `GencanPack.with_restart` at 2.0 Å.

`pack_peo_mesh.py` runs `LatticeGrow` alone. Its cavity is the shipped
`examples/pack_peo/dendrite.stl` — a watertight dendrite, a trunk that forks
three times into 29 branches, 8 732 triangles, authored for a 130 Å cell with a
≈362 000 Å³ cavity — loaded with `scale` so one mesh serves any cell size. The
default 200 × EO25 (35 600 atoms) fills it at 1.026 g/cm³ in 6 s, with 3% of
atoms — hydrogens and end groups on backbone atoms next to the wall — up to
1.3 Å outside the surface.

Three things that scene makes concrete:

- **The push-off is the wrong follow-up at melt density in a cavity.** It has
  nowhere to put the overlap it resolves except through the wall: on this scene
  `GencanPack.with_restart` spent 1 h 45 min moving `fdist` 3.99 → 3.28 while
  `frest` went 0.42 → 6.43 — worst excursion 6.5 Å → 25 Å. Residual contacts at
  melt density are honest, and the force field downstream removes them.
- **A template owes the grower its topology, not its geometry.** Bond lengths,
  angles and torsions are one conformer of that topology, and the force field
  downstream sets them in its first steps. So the backbone atoms *are* the
  lattice sites: the mask's confinement and the guard's self-avoidance are the
  molecule's, torsions are exactly the trans/gauche± the prior drew, angles are
  the lattice's 109.471°, and bond lengths are the lattice step, which is sized
  from the template's own mean backbone bond. Only hydrogens and side atoms
  hang off with template geometry, so a mesh that has to hold a wall is
  authored with about a bond length of clearance in it.
- **`precision` is a distance in disguise.** `frest` is the largest per-atom
  `0.01 · d²`, so `frest < precision` means `d < 10·√precision`: the default
  `1e-2` calls a run converged with an atom 1 Å outside a region.

Each example writes its outputs to `python/examples/out/` (created on
demand, git-ignored) — the path is script-relative, so the working
directory does not matter:

- `{stem}.mrec` — molrs scientific record (`molrs.io.write_mrec_frame`)
- `{stem}.lammpstrj` — LAMMPS dump custom (OVITO particle topology)
- `{stem}.dump.local` — LAMMPS dump local bonds (`batom1`/`batom2`), for
  OVITO [Load trajectory](https://www.ovito.org/manual/reference/pipelines/modifiers/load_trajectory.html)
  (open the `.lammpstrj`, then overlay the dump local file). Skipped when
  the packed frame has no bonds (the water-cube template).

A molrs region answers the two region questions for a batch of points, so
a caller can check what the packer was told to enforce:

```python
cavity = molrs.core.Polyhedron(molrs.io.read_stl("dendrite.stl"))
inside = cavity.contains(state.positions)   # (n,) bool
depth = cavity.distance(state.positions)    # (n,) Å, negative inside
```

`pack_peo_void.py` needs no mesh at all: it reads a LAMMPS data file with
`molrs.io.read_lammps_data`, takes the bonded atoms as the polymer, builds
`molrs.core.SphereUnion(centers, bead_radius + probe, box=...)` and grows PEO
inside `~polymer` — the solvent-accessible void, periodic like the frame.

## Example: mixture

The `pack_mixture.py` example reproduces Packmol's classic `mixture.inp`:

```python
import molrs
from molpack import GencanPack, Target

water_frame = molrs.io.read_pdb("water.pdb")
urea_frame  = molrs.io.read_pdb("urea.pdb")

box = molrs.core.Cuboid([0, 0, 0], [40, 40, 40])

water = Target(water_frame, count=1000).with_name("water").with_restraint(box)
urea  = Target(urea_frame,  count=400).with_name("urea").with_restraint(box)

packer = GencanPack().with_tolerance(2.0).with_seed(1_234_567)
result = packer.run([water, urea], max_loops=400)
print(f"converged={result.converged}  natoms={result.natoms}")
```

## Example: water cube

```python
import molrs
import numpy as np
from molpack import GencanPack, Target

frame = molrs.core.Frame({
    "atoms": {
        "x": np.array([0.00,  0.9572, -0.2400]),
        "y": np.array([0.00,  0.0000,  0.9266]),
        "z": np.zeros(3),
        "element": ["O", "H", "H"],
    }
})

water = Target(frame, count=100).with_name("water").with_restraint(
    molrs.core.Cuboid([0, 0, 0], [30, 30, 30])
)
packer = GencanPack().with_tolerance(2.0).with_progress(False).with_seed(42)
result = packer.run([water], max_loops=200)
print(f"converged={result.converged}  natoms={result.natoms}")
```
