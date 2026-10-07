"""Packmol mixture example: water + urea co-packed in a 40 Å box.

Equivalent to Packmol's ``mixture.inp`` and the Rust
``pack_mixture`` example.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs
import numpy as np

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_mixture"
OUT = HERE / "out"


def main() -> None:
    water_frame = molrs.io.read_pdb(str(DATA / "water.pdb"))
    urea_frame = molrs.io.read_pdb(str(DATA / "urea.pdb"))

    box = molrs.core.Cuboid([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])

    water = (
        molpack.Target(water_frame, count=1000).with_name("water").with_restraint(box)
    )
    urea = molpack.Target(urea_frame, count=400).with_name("urea").with_restraint(box)

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.GencanPack().with_progress(show_progress)

    result = packer.run([water, urea], max_loops=400)

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
    molrs.io.write_mrec_frame(str(OUT / "pack_mixture.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_mixture.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(str(OUT / "pack_mixture.dump.local"), [packed])


if __name__ == "__main__":
    main()
