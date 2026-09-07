"""Linear PEO grown inside a branched STL cavity.

``examples/pack_peo/dendrite.stl`` is a watertight dendrite — a trunk
forking three times into 29 branches, 8732 triangles, authored for a
130 Å cell whose cavity is ≈362 000 Å³. The default 200 × EO25 fills that
at ≈1.03 g/cm³, PEO melt density.
``StlRegion.from_file`` loads it (``scale`` maps the file to whatever
``edge`` you ask for), ``Target.with_restraint`` confines the chains to
it, and ``LatticeGrow`` at 2.0 Å walks Region ∩ lattice — diamond sites
outside the mesh are blocked.

Growth is the whole pipeline here. A seeded push-off
(``GenCanPack().with_restart(grown)``) is the right follow-up for a dilute
box, but at melt density in a cavity it has nowhere to put the overlap it
resolves except through the wall: on this scene it spent 1 h 45 min to
move ``fdist`` 3.99 → 3.28 while ``frest`` went 0.42 → 6.43, i.e. the
worst excursion grew from 6.5 Å to 25 Å. Residual contacts at melt
density are honest, and the force field downstream is what removes them.

A template owes the grower its topology, not its geometry: the backbone
atoms *are* the masked lattice sites, so confinement is the molecule's and
not merely the walk's. Only hydrogens and side atoms hang off with the
template's local geometry, so a mesh that has to hold a wall is authored
with about a bond length of clearance in it.

The mesh path is printed at the start and the frame lands in
``python/examples/out/`` — drop both into a viewer to see the chains inside
the cavity. Passing your own mesh keeps ``scale = 1``, so author it for the
``edge`` you intend to pass.

::

    python python/examples/pack_peo_stl.py 25 200 130 42
    python python/examples/pack_peo_stl.py 25 200 130 42 cavity.stl
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import molpy as mp
import molrs
import numpy as np
from molpy.builder.assembly import (
    MonomerLibrary,
    PolymerBuilder,
    ResiduePlacer,
    SiteMap,
)
from molpy.conformer import Conformer
from molpy.core.atomistic import Atomistic

import molpack

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MESH = HERE.parent.parent / "examples" / "pack_peo" / "dendrite.stl"
#: Cell the shipped mesh is authored for, in Å.
MESH_EDGE = 130.0
ETHER = "[O;%a:1][H].[C:2][O;%b][H]>>[O:1][C:2]"
PEO_C_INF = 5.5
TET = 1.9106332


def make_linear(n: int, *, seed: int = 42) -> Atomistic:
    eo, _ = Conformer(add_hydrogens=True, seed=seed).generate(mp.io.read_smiles("OCCO"))
    SiteMap(eo).label_elements("O", "a", "b")
    return PolymerBuilder(
        MonomerLibrary({"EO": eo}),
        mp.Reaction(ETHER),
        placer=ResiduePlacer(),
    ).build_linear("EO", n)


def pack_stl(
    n: int,
    n_mol: int,
    edge: float,
    seed: int,
    stl_path: Path = MESH,
    scale: float = 1.0,
) -> molpack.State:
    frame = make_linear(n, seed=seed).to_frame()
    cavity = molpack.StlRegion.from_file(stl_path, scale)
    target = (
        molpack.Target(frame, n_mol)
        .with_name("lin-PEO")
        .with_restraint(cavity)
        .with_atom_radius(
            [
                i
                for i, e in enumerate(frame["atoms"].view("element"))
                if str(e).strip() == "H"
            ],
            0.2,
        )
    )
    prior = molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)
    grown = (
        molpack.LatticeGrow(prior)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_periodic_box([0.0, 0.0, 0.0], [edge, edge, edge])
        .run([target], max_loops=max(40, n_mol * 8))
    )
    print(
        f"  grow         : converged={grown.converged}  "
        f"fdist={grown.fdist:.4e}  frest={grown.frest:.4e}  "
        f"degraded={grown.degraded}"
    )
    return grown


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    n = int(args[0]) if args else 25
    n_mol = int(args[1]) if len(args) > 1 else 200
    edge = float(args[2]) if len(args) > 2 else MESH_EDGE
    seed = int(args[3]) if len(args) > 3 else 42
    mesh = Path(args[4]) if len(args) > 4 else MESH
    scale = edge / MESH_EDGE if mesh == MESH else 1.0
    print("── STL-confined linear PEO (LatticeGrow inside StlRegion) ──")
    print(f"  copies       : {n_mol} × EO{n}   cell [0, {edge}]³ Å")
    print(f"  mesh         : {mesh}  (scale {scale:g})")
    t0 = time.perf_counter()
    state = pack_stl(n, n_mol, edge, seed, mesh, scale)
    packed = state.frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.mrec.write_frame(str(OUT / "pack_peo_stl.mrec"), packed)
    molrs.io.write_lammps_traj(str(OUT / "pack_peo_stl.lammpstrj"), [packed])
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(str(OUT / "pack_peo_stl.dump.local"), [packed])
    print(f"  wall         : {time.perf_counter() - t0:.3f} s")


if __name__ == "__main__":
    main()
