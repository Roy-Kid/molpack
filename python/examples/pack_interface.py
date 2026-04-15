"""Packmol interface example: water / chloroform interface + a fixed molecule.

Reproduces Packmol's ``interface.inp``. The ``t3`` molecule is pinned
at the origin with explicit Euler angles — fixed placement.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpack
from _common import read_pdb_as_arrays

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_interface"


def main() -> None:
    water_pos, water_rad, water_els = read_pdb_as_arrays(DATA / "water.pdb")
    chlor_pos, chlor_rad, chlor_els = read_pdb_as_arrays(DATA / "chloroform.pdb")
    t3_pos, t3_rad, t3_els = read_pdb_as_arrays(DATA / "t3.pdb")

    water = (
        molpack.Target.from_coords(water_pos, water_rad, count=100, elements=water_els)
        .with_name("water")
        .with_constraint(molpack.InsideBox([-20.0, 0.0, 0.0], [0.0, 39.0, 39.0]))
    )
    chloroform = (
        molpack.Target.from_coords(chlor_pos, chlor_rad, count=30, elements=chlor_els)
        .with_name("chloroform")
        .with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [21.0, 39.0, 39.0]))
    )
    t3 = (
        molpack.Target.from_coords(t3_pos, t3_rad, count=1, elements=t3_els)
        .with_name("t3")
        .with_center()
        .fixed_at_with_euler(
            position=[0.0, 20.0, 20.0],
            euler=[1.57, 1.57, 1.57],
        )
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack([water, chloroform, t3], max_loops=400, seed=1_234_567)

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
