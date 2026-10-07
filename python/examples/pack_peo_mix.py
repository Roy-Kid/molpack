"""Open-space mixed PEO: linear chains and 4-arm stars in one box.

Two ``Target``s, one ``LatticeGrow.run``, then ``GencanPack.with_restart``.
The box is density-sized from the total mass of both species. Rings stay
out of this scene — both growers raise ``RingTemplate``.

::

    python python/examples/pack_peo_mix.py 4 2 4 4 0.5 42
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import molrs
import numpy as np
from molrs.core import Atomistic

import molpack

OUT = Path(__file__).resolve().parent / "out"

PEO_C_INF = 5.5
TET = 1.9106332
N_ARMS = 4


EO_UNIT = "[<]OCC[>]"  # -O-CH2-CH2-, ports on O (<) and C (>)
CORE_UNIT = "C(C[>])(C[>])(C[>])C[>]"  # pentaerythritol-like four-arm core


def _unit(name: str, body: str, seed: int) -> Atomistic:
    """One CGsmiles unit with its ports, as a 3D molecule with hydrogens."""
    template = molrs.io.smiles.SmilesIr.from_fragment(body).to_template()
    return molrs.conformer.Conformer(seed=seed).generate(template)[0]


def _grow(topology: str, library: dict[str, Atomistic]) -> Atomistic:
    """Grow the CGsmiles ``topology`` from ``library`` into one molecule."""
    sites = molrs.io.cgsmiles.CgSmilesIr(topology).to_coarsegrain()
    return molrs.builder.Assembler(library, molrs.builder.GrowthPlacer()).assemble(
        sites, Atomistic
    )


def linear_topology(n: int) -> str:
    return f"{{[#EO]|{n}}}"


def make_linear(n: int, *, seed: int = 42) -> Atomistic:
    if n < 1:
        raise ValueError(f"linear n must be >= 1, got {n}")
    return _grow(linear_topology(n), {"EO": _unit("EO", EO_UNIT, seed)})


def star_topology(arm_length: int) -> str:
    arm = "[#EO]" * arm_length
    return "{[#X4]" + f"({arm})" * (N_ARMS - 1) + arm + "}"


def make_star(arm_length: int, *, seed: int = 42) -> Atomistic:
    if arm_length < 1:
        raise ValueError(f"star arm_length must be >= 1, got {arm_length}")
    library = {"EO": _unit("EO", EO_UNIT, seed), "X4": _unit("X4", CORE_UNIT, seed + 1)}
    return _grow(star_topology(arm_length), library)


def _report_graph(polymer: Atomistic) -> None:
    n_at = polymer.n_atoms
    n_bd = len(list(polymer.bonds))
    shape = "tree" if n_bd == n_at - 1 else ("unicyclic" if n_bd == n_at else "cyclic+")
    print(f"  template     : {n_at} atoms, {n_bd} bonds  ({shape})")


def _target(polymer: Atomistic, n_mol: int, name: str) -> molpack.Target:
    frame = polymer.to_frame()
    target = molpack.Target(frame, n_mol).with_name(name)
    raw = os.environ.get("PEO_H_RADIUS", "0.2")
    if raw not in ("", "off", "none"):
        elems = list(frame["atoms"]["element"])
        h_idx = [i for i, e in enumerate(elems) if str(e).strip() == "H"]
        if h_idx:
            r = float(raw)
            target = target.with_atom_radius(h_idx, r)
            print(f"  H radius     : {r} Å  ({len(h_idx)} hydrogens)")
    return target


def pack_mix(
    n: int,
    arm_length: int,
    n_linear: int,
    n_star: int,
    density: float,
    seed: int,
):
    if n_linear < 1 or n_star < 1:
        raise ValueError("mix needs at least one linear copy and one star copy")
    linear = make_linear(n, seed=seed)
    star = make_star(arm_length, seed=seed)
    print(f"  linear       : {linear_topology(n)}  × {n_linear}")
    _report_graph(linear)
    print(f"  star         : {star_topology(arm_length)}  × {n_star}")
    _report_graph(star)
    print(
        f"  copies       : {n_linear} linear + {n_star} star   density {density} g/cm³"
    )
    targets = [
        _target(linear, n_linear, "lin-PEO"),
        _target(star, n_star, "star-PEO"),
    ]
    progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "0") != "0"
    prior = molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)
    print(
        "  lattice      : LatticeGrow occupancy-guard @ 2.0 Å → "
        "GencanPack.with_restart @ 2.0 Å"
    )
    grown = (
        molpack.LatticeGrow(prior)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_density(density)
        .with_progress(progress)
        .run(targets, max_loops=max(40, (n_linear + n_star) * 8))
    )
    print(
        f"  grow         : converged={grown.converged}  "
        f"fdist={grown.fdist:.4e}  degraded={grown.degraded}  "
        f"intra scored {grown.intra.scored:.3f} Å"
    )
    pushed = (
        molpack.GencanPack()
        .with_restart(grown)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_progress(progress)
        .run(targets, max_loops=80)
    )
    print(f"  push-off     : converged={pushed.converged}  fdist={pushed.fdist:.4e}")
    return pushed


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    n = int(args[0]) if args else 4
    arm_length = int(args[1]) if len(args) > 1 else 2
    n_linear = int(args[2]) if len(args) > 2 else 4
    n_star = int(args[3]) if len(args) > 3 else 4
    density = float(args[4]) if len(args) > 4 else 0.5
    seed = int(args[5]) if len(args) > 5 else 42
    print("── mixed PEO (linear + 4-arm star, one box) ──")
    t0 = time.perf_counter()
    packed = pack_mix(n, arm_length, n_linear, n_star, density, seed).frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.core.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.write_mrec_frame(str(OUT / "pack_peo_mix.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_peo_mix.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(str(OUT / "pack_peo_mix.dump.local"), [packed])
    print(f"  wall         : {time.perf_counter() - t0:.3f} s")


if __name__ == "__main__":
    main()
