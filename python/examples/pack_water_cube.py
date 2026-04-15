"""Pack 100 water molecules into a 30x30x30 cubic box.

Minimal standalone example — no PDB file required.
"""

from __future__ import annotations

import numpy as np

import molpack


def main() -> None:
    # Water geometry: O at origin, two Hs 0.96 Å away.
    frame = {
        "atoms": {
            "x": np.array([0.0, 0.9572, -0.2400], dtype=np.float64),
            "y": np.array([0.0, 0.0, 0.9266], dtype=np.float64),
            "z": np.zeros(3, dtype=np.float64),
            "element": ["O", "H", "H"],
        }
    }

    water = molpack.Target("water", frame, count=100).with_restraint(
        molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
    )

    packer = molpack.Molpack(tolerance=2.0, precision=0.01, progress=False)
    result = packer.pack([water], max_loops=200, seed=42)

    print(f"converged = {result.converged}")
    print(f"natoms    = {result.natoms}")
    print(f"fdist     = {result.fdist:.4f}")
    print(f"frest     = {result.frest:.4f}")
    print(f"positions shape = {result.positions.shape}")


if __name__ == "__main__":
    main()
