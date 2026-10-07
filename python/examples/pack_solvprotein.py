"""Packmol solvprotein example: one fixed solute + water + ions in a 50 Å sphere.

Reproduces Packmol's ``solvprotein.inp``. The solute is centred at
the origin and pinned; water and monoatomic ions are packed into a
spherical shell around it.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs
import numpy as np

import molpack
from molpack import CenteringMode

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_solvprotein"
OUT = HERE / "out"


def main() -> None:
    protein_frame = molrs.io.read_pdb(str(DATA / "protein.pdb"))
    water_frame = molrs.io.read_pdb(str(DATA / "water.pdb"))
    sodium_frame = molrs.io.read_pdb(str(DATA / "sodium.pdb"))
    chloride_frame = molrs.io.read_pdb(str(DATA / "chloride.pdb"))

    sphere = molrs.core.Sphere([0.0, 0.0, 0.0], 50.0)

    protein = (
        molpack.Target(protein_frame, count=1)
        .with_name("protein")
        .with_centering(CenteringMode.CENTER)
        .fixed_at([0.0, 0.0, 0.0])
    )
    water = (
        molpack.Target(water_frame, count=1000)
        .with_name("water")
        .with_restraint(sphere)
    )
    sodium = (
        molpack.Target(sodium_frame, count=30)
        .with_name("sodium")
        .with_restraint(sphere)
    )
    chloride = (
        molpack.Target(chloride_frame, count=20)
        .with_name("chloride")
        .with_restraint(sphere)
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.GenCanPack().with_progress(show_progress)

    result = packer.run(
        [protein, water, sodium, chloride],
        max_loops=800,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )
    packed = result.frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.core.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.write_mrec_frame(str(OUT / "pack_solvprotein.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_solvprotein.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(
            str(OUT / "pack_solvprotein.dump.local"), [packed]
        )


if __name__ == "__main__":
    main()
