"""Packmol solvprotein example: fixed protein + water + ions in a 50 Å sphere.

Reproduces Packmol's ``solvprotein.inp``. The protein is centred at
the origin and pinned; water and monoatomic ions are packed into a
spherical shell around it.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpack
from _common import read_pdb_as_arrays

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_solvprotein"


def main() -> None:
    protein_pos, protein_rad, protein_els = read_pdb_as_arrays(DATA / "protein.pdb")
    water_pos, water_rad, water_els = read_pdb_as_arrays(DATA / "water.pdb")
    sodium_pos, sodium_rad, sodium_els = read_pdb_as_arrays(DATA / "sodium.pdb")
    chlor_pos, chlor_rad, chlor_els = read_pdb_as_arrays(DATA / "chloride.pdb")

    sphere = molpack.InsideSphere(50.0, [0.0, 0.0, 0.0])

    protein = (
        molpack.Target.from_coords(
            protein_pos, protein_rad, count=1, elements=protein_els
        )
        .with_name("protein")
        .with_center()
        .fixed_at([0.0, 0.0, 0.0])
    )
    water = (
        molpack.Target.from_coords(water_pos, water_rad, count=1000, elements=water_els)
        .with_name("water")
        .with_constraint(sphere)
    )
    sodium = (
        molpack.Target.from_coords(
            sodium_pos, sodium_rad, count=30, elements=sodium_els
        )
        .with_name("sodium")
        .with_constraint(sphere)
    )
    chloride = (
        molpack.Target.from_coords(chlor_pos, chlor_rad, count=20, elements=chlor_els)
        .with_name("chloride")
        .with_constraint(sphere)
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(show_progress)

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
