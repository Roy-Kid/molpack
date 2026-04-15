"""Packmol solvprotein example: fixed protein + water + ions in a 50 Å sphere.

Reproduces Packmol's ``solvprotein.inp``. The protein is centred at
the origin and pinned; water and monoatomic ions are packed into a
spherical shell around it.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_solvprotein"


def main() -> None:
    protein_frame = molrs.read_pdb(str(DATA / "protein.pdb"))
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    sodium_frame = molrs.read_pdb(str(DATA / "sodium.pdb"))
    chloride_frame = molrs.read_pdb(str(DATA / "chloride.pdb"))

    sphere = molpack.InsideSphere(50.0, [0.0, 0.0, 0.0])

    protein = (
        molpack.Target("protein", protein_frame, count=1)
        .with_center()
        .fixed_at([0.0, 0.0, 0.0])
    )
    water = molpack.Target("water", water_frame, count=1000).with_restraint(sphere)
    sodium = molpack.Target("sodium", sodium_frame, count=30).with_restraint(sphere)
    chloride = molpack.Target("chloride", chloride_frame, count=20).with_restraint(
        sphere
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Molpack(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack(
        [protein, water, sodium, chloride],
        max_loops=800,
        seed=1_234_567,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
