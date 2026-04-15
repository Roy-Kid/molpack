"""Packmol mixture example: water + urea co-packed in a 40 Å box.

Equivalent to Packmol's ``mixture.inp`` and the Rust
``pack_mixture`` example.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpack
from _common import read_pdb_as_arrays

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_mixture"


def main() -> None:
    water_pos, water_rad, water_els = read_pdb_as_arrays(DATA / "water.pdb")
    urea_pos, urea_rad, urea_els = read_pdb_as_arrays(DATA / "urea.pdb")

    box = molpack.InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])

    water = (
        molpack.Target.from_coords(water_pos, water_rad, count=1000, elements=water_els)
        .with_name("water")
        .with_constraint(box)
    )
    urea = (
        molpack.Target.from_coords(urea_pos, urea_rad, count=400, elements=urea_els)
        .with_name("urea")
        .with_constraint(box)
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack([water, urea], max_loops=400, seed=1_234_567)

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
