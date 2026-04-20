"""Packmol interface example: water / chloroform interface + a fixed molecule.

Reproduces Packmol's ``interface.inp``. The ``t3`` molecule is pinned
at the origin with explicit Euler angles — fixed placement.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack
from molpack import Angle, CenteringMode

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_interface"


def main() -> None:
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    chlor_frame = molrs.read_pdb(str(DATA / "chloroform.pdb"))
    t3_frame = molrs.read_pdb(str(DATA / "t3.pdb"))

    water = (
        molpack.Target(water_frame, count=100)
        .with_name("water")
        .with_restraint(
            molpack.InsideBoxRestraint([-20.0, 0.0, 0.0], [0.0, 39.0, 39.0])
        )
    )
    chloroform = (
        molpack.Target(chlor_frame, count=30)
        .with_name("chloroform")
        .with_restraint(molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [21.0, 39.0, 39.0]))
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

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = (
        molpack.Molpack()
        .with_tolerance(2.0)
        .with_precision(0.01)
        .with_progress(show_progress)
        .with_seed(1_234_567)
    )

    result = packer.pack([water, chloroform, t3], max_loops=400)

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
