"""Linear PEO grown inside a branched STL cavity.

``examples/pack_peo/dendrite.stl`` is a watertight dendrite — a trunk
forking three times into 29 branches, 8732 triangles, authored for a
130 Å cell whose cavity is ≈362 000 Å³. The default 200 × EO25 fills that
at ≈1.03 g/cm³, PEO melt density.
``StlRegion.from_file`` loads it (``scale`` maps the file to whatever
``edge`` you ask for), ``Target.with_restraint`` confines the chains to
it, and ``LatticeGrow`` at 2.0 Å walks Region ∩ lattice (diamond sites
outside the mesh are blocked) before ``GenCanPack.with_restart`` pushes
the contacts open at ``precision=1e-4`` — the default 1e-2 lets an atom
sit ~1 Å outside and still call the run converged.

The mask confines the *walk*, not the decorated atoms: decoration
rebuilds the molecule from its own bonds and angles along that track and
drifts off it, by more the longer the backbone. ``with_track_tweak``
bounds how far a torsion may bend to follow the track and is the lever on
that drift — at EO25 the worst excursion after the walk is ≈36 Å at 0.0,
≈16 Å at the 0.35 default, ≈3 Å at 1.5, and flat beyond. A mesh that has
to hold a wall is still authored with the clearance already in it.

Drop the mesh into a viewer next to the packed frame to see the cavity:
both paths are printed at the end.

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
            0.85,
        )
    )
    prior = molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)
    grown = (
        molpack.LatticeGrow(prior)
        # A 25-mer backbone is 75 atoms long, and decoration rebuilds it from
        # the template's own bonds and angles along the walk: at the 0.35 rad
        # default the error compounds to ~16 Å by the last residue. 1.5 rad
        # keeps it near 3 Å. The price is that a hooked torsion can then leave
        # its RIS state entirely, so the realized statistics are the lattice's
        # rather than the prior's — drop back to the default when the
        # conformer matters more than the wall.
        .with_track_tweak(1.5)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_periodic_box([0.0, 0.0, 0.0], [edge, edge, edge])
        .run([target], max_loops=max(40, n_mol * 8))
    )
    print(
        f"  grow         : converged={grown.converged}  "
        f"fdist={grown.fdist:.4e}  frest={grown.frest:.4e}  "
        f"softened={grown.softened}"
    )
    pushed = (
        molpack.GenCanPack()
        .with_restart(grown)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_precision(1e-4)
        .run([target], max_loops=80)
    )
    print(
        f"  push-off     : converged={pushed.converged}  "
        f"fdist={pushed.fdist:.4e}  frest={pushed.frest:.4e}"
    )
    return pushed


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
