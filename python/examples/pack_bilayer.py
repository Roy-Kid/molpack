"""Packmol bilayer example: a double layer with water above and below.

Based on Packmol's ``bilayer.inp`` — atom-subset restraints pin the
head atoms and the tail atoms of each molecule into their
respective slabs.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs
import numpy as np

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_bilayer"
OUT = HERE / "out"


def main() -> None:
    water_frame = molrs.io.read_pdb(str(DATA / "water.pdb"))
    lipid_frame = molrs.io.read_pdb(str(DATA / "palmitoil.pdb"))

    water_low = (
        molpack.Target(water_frame, count=50)
        .with_name("water_low")
        .with_restraint(molrs.core.Cuboid([0.0, 0.0, -10.0], [40.0, 40.0, 10.0]))
    )

    water_high = (
        molpack.Target(water_frame, count=50)
        .with_name("water_high")
        .with_restraint(molrs.core.Cuboid([0.0, 0.0, 28.0], [40.0, 40.0, 10.0]))
    )

    lipid_low = (
        molpack.Target(lipid_frame, count=10)
        .with_name("lipid_low")
        .with_restraint(molrs.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 14.0]))
        # 0-based: Packmol .inp atoms 32/33 → indices 31/32 for tails below z=2
        .with_atom_restraint(
            [30, 31], molrs.core.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 2.0])
        )
        # Packmol .inp atoms 1/2 → indices 0/1 for heads above z=12
        .with_atom_restraint(
            [0, 1], ~molrs.core.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 12.0])
        )
    )

    lipid_high = (
        molpack.Target(lipid_frame, count=10)
        .with_name("lipid_high")
        .with_restraint(molrs.core.Cuboid([0.0, 0.0, 14.0], [40.0, 40.0, 14.0]))
        # heads below z=16
        .with_atom_restraint(
            [0, 1], molrs.core.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 16.0])
        )
        # tails above z=26
        .with_atom_restraint(
            [30, 31], ~molrs.core.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 26.0])
        )
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.GencanPack().with_progress(show_progress)

    result = packer.run(
        [water_low, water_high, lipid_low, lipid_high],
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
    molrs.io.write_mrec_frame(str(OUT / "pack_bilayer.mrec"), packed)
    molrs.io.write_lammps_dump_trajectory(
        str(OUT / "pack_bilayer.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].n_rows:
        molrs.io.write_lammps_dump_local(str(OUT / "pack_bilayer.dump.local"), [packed])


if __name__ == "__main__":
    main()
