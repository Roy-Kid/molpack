"""Pack 100 water molecules into a 30x30x30 cubic box.

Minimal example — builds the template with ``molrs.store.Frame`` (no PDB
file), so it needs ``molcrafts-molrs`` but no structure files on disk.
"""

from __future__ import annotations

from pathlib import Path

import molrs
import numpy as np

import molpack

OUT = Path(__file__).resolve().parent / "out"


def main() -> None:
    # Water geometry: O at origin, two Hs 0.96 Å away.
    frame = molrs.store.Frame(
        {
            "atoms": {
                "x": np.array([0.0, 0.9572, -0.2400], dtype=np.float64),
                "y": np.array([0.0, 0.0, 0.9266], dtype=np.float64),
                "z": np.zeros(3, dtype=np.float64),
                "element": ["O", "H", "H"],
            }
        }
    )

    water = (
        molpack.Target(frame, count=100)
        .with_name("water")
        .with_restraint(molrs.spatial.Cuboid([0.0, 0.0, 0.0], [30.0, 30.0, 30.0]))
    )

    packer = molpack.GenCanPack()
    result = packer.run([water], max_loops=200)

    print(f"converged = {result.converged}")
    print(f"natoms    = {result.natoms}")
    print(f"fdist     = {result.fdist:.4f}")
    print(f"frest     = {result.frest:.4f}")
    print(f"positions shape = {result.positions.shape}")
    packed = result.frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.spatial.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.write_mrec(str(OUT / "pack_water_cube.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_water_cube.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(
            str(OUT / "pack_water_cube.dump.local"), [packed]
        )


if __name__ == "__main__":
    main()
