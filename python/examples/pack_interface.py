"""Packmol interface example: water / chloroform interface + a fixed molecule.

Reproduces Packmol's ``interface.inp``. The ``t3`` molecule is pinned
at the origin with explicit Euler angles — fixed placement.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_interface"


def main() -> None:
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    chlor_frame = molrs.read_pdb(str(DATA / "chloroform.pdb"))
    t3_frame = molrs.read_pdb(str(DATA / "t3.pdb"))

    water = molpack.Target("water", water_frame, count=100).with_restraint(
        molpack.InsideBox([-20.0, 0.0, 0.0], [0.0, 39.0, 39.0])
    )
    chloroform = molpack.Target("chloroform", chlor_frame, count=30).with_restraint(
        molpack.InsideBox([0.0, 0.0, 0.0], [21.0, 39.0, 39.0])
    )
    t3 = (
        molpack.Target("t3", t3_frame, count=1)
        .with_center()
        .fixed_at_with_euler(
            position=[0.0, 20.0, 20.0],
            euler=[1.57, 1.57, 1.57],
        )
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Molpack(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack([water, chloroform, t3], max_loops=400, seed=1_234_567)

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
