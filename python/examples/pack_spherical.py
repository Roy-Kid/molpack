"""Packmol spherical example: nested double-layered shell.

Reproduces Packmol's ``spherical.inp``: four shells around the origin,
with atom-subset restraints pinning molecule heads and tails onto the
correct inner/outer surfaces of the double layer.
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs
import numpy as np

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_spherical"
OUT = HERE / "out"

ORIGIN = [0.0, 0.0, 0.0]


def main() -> None:
    water_frame = molrs.io.read_pdb(str(DATA / "water.pdb"))
    lipid_frame = molrs.io.read_pdb(str(DATA / "palmitoil.pdb"))

    # 1. Inner water sphere (r = 13).
    water_inner = (
        molpack.Target(water_frame, count=308)
        .with_name("water_inner")
        .with_restraint(molrs.core.Sphere(ORIGIN, 13.0))
    )

    # 2. Inner layer: head atom (0-based index 36) inside r=14,
    # tail atom (0-based index 4) outside r=26.
    lipid_inner = (
        molpack.Target(lipid_frame, count=90)
        .with_name("lipid_inner")
        .with_atom_restraint([36], molrs.core.Sphere(ORIGIN, 14.0))
        .with_atom_restraint([4], ~molrs.core.Sphere(ORIGIN, 26.0))
    )

    # 3. Outer layer: tail atom 4 inside r=29, head atom 36 outside r=41.
    lipid_outer = (
        molpack.Target(lipid_frame, count=300)
        .with_name("lipid_outer")
        .with_atom_restraint([4], molrs.core.Sphere(ORIGIN, 29.0))
        .with_atom_restraint([36], ~molrs.core.Sphere(ORIGIN, 41.0))
    )

    # 4. Outer water shell: inside ±47.5 box and outside sphere r=43.
    water_outer = (
        molpack.Target(water_frame, count=17_536)
        .with_name("water_outer")
        .with_restraint(molrs.core.Cuboid([-47.5, -47.5, -47.5], [95.0, 95.0, 95.0]))
        .with_restraint(~molrs.core.Sphere(ORIGIN, 43.0))
    )

    log_level = (
        "progress"
        if os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
        else "quiet"
    )
    packer = molpack.GencanPack().with_log_level(log_level)

    result = packer.run(
        [water_inner, lipid_inner, lipid_outer, water_outer],
        max_loops=800,
    )

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
    molrs.io.write_mrec_frame(str(OUT / "pack_spherical.mrec"), packed)
    molrs.io.write_lammps_dump_trajectory(
        str(OUT / "pack_spherical.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].n_rows:
        molrs.io.write_lammps_dump_local(
            str(OUT / "pack_spherical.dump.local"), [packed]
        )


if __name__ == "__main__":
    main()
