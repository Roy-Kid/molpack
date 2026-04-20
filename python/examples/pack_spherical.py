"""Packmol spherical example: nested double-layered vesicle.

Reproduces Packmol's ``spherical.inp``: four shells around the origin,
with atom-subset restraints pinning lipid heads and tails onto the
correct inner/outer surfaces of the bilayer.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_spherical"

ORIGIN = [0.0, 0.0, 0.0]


def main() -> None:
    water_frame = molrs.read_pdb(str(DATA / "water.pdb"))
    lipid_frame = molrs.read_pdb(str(DATA / "palmitoil.pdb"))

    # 1. Inner water sphere (r = 13).
    water_inner = (
        molpack.Target(water_frame, count=308)
        .with_name("water_inner")
        .with_restraint(molpack.InsideSphereRestraint(ORIGIN, 13.0))
    )

    # 2. Inner lipid layer: head atom (0-based index 36) inside r=14,
    # tail atom (0-based index 4) outside r=26.
    lipid_inner = (
        molpack.Target(lipid_frame, count=90)
        .with_name("lipid_inner")
        .with_atom_restraint([36], molpack.InsideSphereRestraint(ORIGIN, 14.0))
        .with_atom_restraint([4], molpack.OutsideSphereRestraint(ORIGIN, 26.0))
    )

    # 3. Outer lipid layer: tail atom 4 inside r=29, head atom 36 outside r=41.
    lipid_outer = (
        molpack.Target(lipid_frame, count=300)
        .with_name("lipid_outer")
        .with_atom_restraint([4], molpack.InsideSphereRestraint(ORIGIN, 29.0))
        .with_atom_restraint([36], molpack.OutsideSphereRestraint(ORIGIN, 41.0))
    )

    # 4. Outer water shell: inside ±47.5 box and outside sphere r=43.
    water_outer = (
        molpack.Target(water_frame, count=17_536)
        .with_name("water_outer")
        .with_restraint(
            molpack.InsideBoxRestraint([-47.5, -47.5, -47.5], [47.5, 47.5, 47.5])
        )
        .with_restraint(molpack.OutsideSphereRestraint(ORIGIN, 43.0))
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
        [water_inner, lipid_inner, lipid_outer, water_outer],
        max_loops=800,
    )

    print(
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )


if __name__ == "__main__":
    main()
