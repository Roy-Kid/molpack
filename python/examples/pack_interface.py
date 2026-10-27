"""Packmol interface example: water / chloroform interface + a fixed molecule.

Reproduces Packmol's ``interface.inp``. The ``t3`` molecule is pinned
at the origin with explicit Euler angles — fixed placement.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpy
import numpy as np

import molpack
from molpack import Angle, CenteringMode

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_interface"
OUT = HERE / "out"


def main() -> None:
    water_frame = molpy.io.read_pdb(DATA / "water.pdb")
    chlor_frame = molpy.io.read_pdb(DATA / "chloroform.pdb")
    t3_frame = molpy.io.read_pdb(DATA / "t3.pdb")

    water = (
        molpack.Target(water_frame, count=100)
        .with_name("water")
        .with_restraint(molpy.core.Cuboid([-20.0, 0.0, 0.0], [20.0, 39.0, 39.0]))
    )
    chloroform = (
        molpack.Target(chlor_frame, count=30)
        .with_name("chloroform")
        .with_restraint(molpy.core.Cuboid([0.0, 0.0, 0.0], [21.0, 39.0, 39.0]))
    )
    t3 = (
        molpack.Target(t3_frame, count=1)
        .with_name("t3")
        .with_centering(CenteringMode.CENTER)
        .fixed_at([0.0, 20.0, 20.0])
        .with_orientation(
            (
                Angle.from_radians(1.57),
                Angle.from_radians(1.57),
                Angle.from_radians(1.57),
            )
        )
    )

    log_level = (
        "progress"
        if os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
        else "quiet"
    )
    packer = molpack.GencanPack().with_log_level(log_level)

    result = packer.run([water, chloroform, t3], max_loops=400)

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
    molpy.io.write_mrec(OUT / "pack_interface.mrec", frame=packed)
    molpy.io.write_lammps_dump(
        OUT / "pack_interface.lammpstrj",
        packed,
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].n_rows:
        molpy.io.write_lammps_dump_local(OUT / "pack_interface.dump.local", [packed])


if __name__ == "__main__":
    main()
