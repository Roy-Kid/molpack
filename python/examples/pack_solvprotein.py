"""Packmol solvprotein example: one fixed solute + water + ions in a 50 Å sphere.

Reproduces Packmol's ``solvprotein.inp``. The solute is centred at
the origin and pinned; water and monoatomic ions are packed into a
spherical shell around it.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpy
import numpy as np

import molpack
from molpack import CenteringMode

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_solvprotein"
OUT = HERE / "out"


def main() -> None:
    protein_frame = molpy.io.read_pdb(DATA / "protein.pdb")
    water_frame = molpy.io.read_pdb(DATA / "water.pdb")
    sodium_frame = molpy.io.read_pdb(DATA / "sodium.pdb")
    chloride_frame = molpy.io.read_pdb(DATA / "chloride.pdb")

    sphere = molpy.core.Sphere([0.0, 0.0, 0.0], 50.0)

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

    log_level = (
        "progress"
        if os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
        else "quiet"
    )
    packer = molpack.GencanPack().with_log_level(log_level)

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
        packed.box = molpy.core.Box.bounding(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            np.ones(3),
            np.array([False, False, False]),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molpy.io.write_mrec(OUT / "pack_solvprotein.mrec", frame=packed)
    molpy.io.write_lammps_dump(
        OUT / "pack_solvprotein.lammpstrj",
        packed,
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].n_rows:
        molpy.io.write_lammps_dump_local(
            str(OUT / "pack_solvprotein.dump.local"), [packed]
        )


if __name__ == "__main__":
    main()
