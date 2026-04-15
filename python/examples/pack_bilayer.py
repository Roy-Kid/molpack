"""Packmol bilayer example: lipid double layer with water above and below.

Based on Packmol's ``bilayer.inp`` — atom-subset restraints pin the
polar heads and the aliphatic tails of each lipid into their
respective slabs.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpack
from _common import read_pdb_as_arrays

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_bilayer"


def main() -> None:
    water_pos, water_rad, water_els = read_pdb_as_arrays(DATA / "water.pdb")
    lipid_pos, lipid_rad, lipid_els = read_pdb_as_arrays(DATA / "palmitoil.pdb")

    water_low = (
        molpack.Target.from_coords(water_pos, water_rad, count=50, elements=water_els)
        .with_name("water_low")
        .with_constraint(molpack.InsideBox([0.0, 0.0, -10.0], [40.0, 40.0, 0.0]))
    )

    water_high = (
        molpack.Target.from_coords(water_pos, water_rad, count=50, elements=water_els)
        .with_name("water_high")
        .with_constraint(molpack.InsideBox([0.0, 0.0, 28.0], [40.0, 40.0, 38.0]))
    )

    lipid_low = (
        molpack.Target.from_coords(lipid_pos, lipid_rad, count=10, elements=lipid_els)
        .with_name("lipid_low")
        .with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 14.0]))
        # tail atoms 31–32 below z=2
        .with_constraint_for_atoms([31, 32], molpack.BelowPlane([0.0, 0.0, 1.0], 2.0))
        # head atoms 1–2 above z=12
        .with_constraint_for_atoms([1, 2], molpack.AbovePlane([0.0, 0.0, 1.0], 12.0))
    )

    lipid_high = (
        molpack.Target.from_coords(lipid_pos, lipid_rad, count=10, elements=lipid_els)
        .with_name("lipid_high")
        .with_constraint(molpack.InsideBox([0.0, 0.0, 14.0], [40.0, 40.0, 28.0]))
        # head atoms 1–2 below z=16
        .with_constraint_for_atoms([1, 2], molpack.BelowPlane([0.0, 0.0, 1.0], 16.0))
        # tail atoms 31–32 above z=26
        .with_constraint_for_atoms([31, 32], molpack.AbovePlane([0.0, 0.0, 1.0], 26.0))
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack(
        [water_low, water_high, lipid_low, lipid_high],
        max_loops=800,
        seed=1_234_567,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
