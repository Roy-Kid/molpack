"""Packmol bilayer example: lipid double layer with water above and below.

Based on Packmol's ``bilayer.inp`` — atom-subset restraints pin the
polar heads and the aliphatic tails of each lipid into their
respective slabs.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_bilayer"


def main() -> None:
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    lipid_frame = molrs.read_pdb(str(DATA / "palmitoil.pdb"))

    water_low = (
        molpack.Target(water_frame, count=50)
        .with_name("water_low")
        .with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, -10.0], [40.0, 40.0, 0.0])
        )
    )

    water_high = (
        molpack.Target(water_frame, count=50)
        .with_name("water_high")
        .with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 28.0], [40.0, 40.0, 38.0])
        )
    )

    lipid_low = (
        molpack.Target(lipid_frame, count=10)
        .with_name("lipid_low")
        .with_restraint(molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 14.0]))
        # 0-based: Packmol .inp atoms 32/33 → indices 31/32 for tails below z=2
        .with_atom_restraint(
            [30, 31], molpack.BelowPlaneRestraint([0.0, 0.0, 1.0], 2.0)
        )
        # Packmol .inp atoms 1/2 → indices 0/1 for heads above z=12
        .with_atom_restraint([0, 1], molpack.AbovePlaneRestraint([0.0, 0.0, 1.0], 12.0))
    )

    lipid_high = (
        molpack.Target(lipid_frame, count=10)
        .with_name("lipid_high")
        .with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 14.0], [40.0, 40.0, 28.0])
        )
        # heads below z=16
        .with_atom_restraint([0, 1], molpack.BelowPlaneRestraint([0.0, 0.0, 1.0], 16.0))
        # tails above z=26
        .with_atom_restraint(
            [30, 31], molpack.AbovePlaneRestraint([0.0, 0.0, 1.0], 26.0)
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

    result = packer.pack(
        [water_low, water_high, lipid_low, lipid_high],
        max_loops=800,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
