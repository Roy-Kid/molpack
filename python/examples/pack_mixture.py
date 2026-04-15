"""Packmol mixture example: water + urea co-packed in a 40 Å box.

Equivalent to Packmol's ``mixture.inp`` and the Rust
``pack_mixture`` example.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_mixture"


def main() -> None:
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    urea_frame = molrs.read_pdb(str(DATA / "urea.pdb"))

    box = molpack.InsideBox([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])

    water = molpack.Target("water", water_frame, count=1000).with_restraint(box)
    urea = molpack.Target("urea", urea_frame, count=400).with_restraint(box)

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Molpack(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack([water, urea], max_loops=400, seed=1_234_567)

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
