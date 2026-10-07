"""Keeping one species from clustering: ions dispersed through a water box.

The packer's pair term stops molecules overlapping and then stops caring, and
nothing in it distinguishes two copies of one species from a copy of each of
two. Ions are therefore free to end up in a corner together as long as they do
not interpenetrate — which is exactly what the first run below shows.

``SelfSeparation`` is the missing statement: *these molecules also keep their
distance from each other*. The second run adds it and nothing else, so the
difference in the printed ion-ion statistics is attributable to that one line.

Run::

    python examples/pack_ion_dispersion.py
"""

from __future__ import annotations

import os
from pathlib import Path

import molrs
import numpy as np

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_solvprotein"
OUT = HERE / "out"

BOX_LO = [0.0, 0.0, 0.0]
BOX_HI = [35.0, 35.0, 35.0]
N_WATER = 400
N_ION = 24
D_MIN = 9.0
SEED = 42


def ion_statistics(pos: np.ndarray) -> dict[str, float]:
    """Ion-ion spacing summary. ``pos`` is (n, 3) of ion positions."""
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    nearest = d.min(axis=1)
    # Mean number of same-species neighbours inside the requested distance —
    # the "am I in a clump" number, 0 when the bound is honoured.
    crowding = float((d < D_MIN).sum(axis=1).mean())
    return {
        "closest pair": float(d.min()),
        "mean nearest neighbour": float(nearest.mean()),
        "neighbours within D_MIN": crowding,
    }


def pack(*, separate: bool, log_level: str) -> tuple[np.ndarray, object]:
    """Pack water + ions once. ``separate`` toggles the one line under test."""
    water_frame = molrs.io.read_pdb(str(DATA / "water.pdb"))
    ion_frame = molrs.io.read_pdb(str(DATA / "sodium.pdb"))
    box = molrs.core.Cuboid(BOX_LO, np.subtract(BOX_HI, BOX_LO))

    water = (
        molpack.Target(water_frame, count=N_WATER)
        .with_name("water")
        .with_restraint(box)
    )
    ions = molpack.Target(ion_frame, count=N_ION).with_name("NA").with_restraint(box)
    if separate:
        ions = ions.with_restraint(molpack.SelfSeparation(D_MIN))

    packer = molpack.GencanPack().with_log_level(log_level).with_seed(SEED)
    result = packer.run([water, ions], max_loops=200)

    # Targets are packed in the order given, so the ions are the trailing
    # `N_ION` atoms (sodium is monatomic).
    pos = np.asarray(result.positions)[-N_ION:]
    return pos, result


def report(label: str, pos: np.ndarray, result) -> None:
    stats = ion_statistics(pos)
    print(f"\n{label}")
    print(
        f"  converged={result.converged} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )
    for name, value in stats.items():
        print(f"  {name:<26} {value:7.2f}")


def write(name: str, result) -> None:
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
    molrs.io.write_mrec_frame(str(OUT / f"{name}.mrec"), packed)
    molrs.io.write_lammps_dump_trajectory(
        str(OUT / f"{name}.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )


def main() -> None:
    log_level = (
        "progress"
        if os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
        else "quiet"
    )
    print(
        f"{N_ION} ions + {N_WATER} water in a "
        f"{BOX_HI[0]:.0f} A box, asking for {D_MIN:.0f} A between ions"
    )

    free_pos, free_result = pack(separate=False, log_level=log_level)
    report("without SelfSeparation (the control)", free_pos, free_result)
    write("pack_ion_dispersion_free", free_result)

    kept_pos, kept_result = pack(separate=True, log_level=log_level)
    report(f"with SelfSeparation({D_MIN})", kept_pos, kept_result)
    write("pack_ion_dispersion_separated", kept_result)

    print(
        "\nThe bound is a lower bound on the ion-ion distance, not a uniformity\n"
        "target: it says how close two ions may come, not how they are spread.\n"
        "Ask for more than fits and the run reports it through frest rather\n"
        "than relaxing it silently."
    )


if __name__ == "__main__":
    main()
