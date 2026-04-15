"""Pack 100 water molecules into a 30x30x30 cubic box.

Minimal standalone example using molpack's Python API.
"""

from __future__ import annotations

import numpy as np

import molpack


def main() -> None:
    # Water geometry: O at origin, two Hs 0.96 Å away.
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.2400, 0.9266, 0.0],
        ],
        dtype=np.float64,
    )
    # VdW radii (Å): O = 1.52, H = 1.20.
    radii = np.array([1.52, 1.20, 1.20], dtype=np.float64)

    water = (
        molpack.Target.from_coords(
            positions, radii, count=100, elements=["O", "H", "H"]
        )
        .with_name("water")
        .with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0]))
    )

    packer = molpack.Packer(tolerance=2.0, precision=0.01).with_progress(False)
    result = packer.pack([water], max_loops=200, seed=42)

    print(f"converged = {result.converged}")
    print(f"natoms    = {result.natoms}")
    print(f"fdist     = {result.fdist:.4f}")
    print(f"frest     = {result.frest:.4f}")
    print(f"positions shape = {result.positions.shape}")


if __name__ == "__main__":
    main()
