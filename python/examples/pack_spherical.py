"""Packmol spherical example: nested double-layered vesicle.

Reproduces Packmol's ``spherical.inp``: four shells around the origin,
with atom-subset restraints pinning lipid heads and tails onto the
correct inner/outer surfaces of the bilayer.
"""

from __future__ import annotations

import os
from pathlib import Path

import molpack
from _common import read_pdb_as_arrays

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_spherical"

ORIGIN = [0.0, 0.0, 0.0]


def main() -> None:
    water_pos, water_rad, water_els = read_pdb_as_arrays(DATA / "water.pdb")
    lipid_pos, lipid_rad, lipid_els = read_pdb_as_arrays(DATA / "palmitoil.pdb")

    # 1. Inner water sphere (r = 13).
    water_inner = (
        molpack.Target.from_coords(water_pos, water_rad, count=308, elements=water_els)
        .with_name("water_inner")
        .with_constraint(molpack.InsideSphere(13.0, ORIGIN))
    )

    # 2. Inner lipid layer: head atom 37 inside r=14, tail atom 5 outside r=26.
    lipid_inner = (
        molpack.Target.from_coords(lipid_pos, lipid_rad, count=90, elements=lipid_els)
        .with_name("lipid_inner")
        .with_constraint_for_atoms([37], molpack.InsideSphere(14.0, ORIGIN))
        .with_constraint_for_atoms([5], molpack.OutsideSphere(26.0, ORIGIN))
    )

    # 3. Outer lipid layer: tail atom 5 inside r=29, head atom 37 outside r=41.
    lipid_outer = (
        molpack.Target.from_coords(lipid_pos, lipid_rad, count=300, elements=lipid_els)
        .with_name("lipid_outer")
        .with_constraint_for_atoms([5], molpack.InsideSphere(29.0, ORIGIN))
        .with_constraint_for_atoms([37], molpack.OutsideSphere(41.0, ORIGIN))
    )

    # 4. Outer water shell: inside ±47.5 box and outside sphere r=43.
    water_outer = (
        molpack.Target.from_coords(
            water_pos, water_rad, count=17_536, elements=water_els
        )
        .with_name("water_outer")
        .with_constraint(
            molpack.InsideBox([-47.5, -47.5, -47.5], [47.5, 47.5, 47.5]).and_(
                molpack.OutsideSphere(43.0, ORIGIN)
            )
        )
    )

    show_progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(show_progress)

    result = packer.pack(
        [water_inner, lipid_inner, lipid_outer, water_outer],
        max_loops=800,
        seed=1_234_567,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
