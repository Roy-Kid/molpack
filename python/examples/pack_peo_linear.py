"""Open-space linear PEO melt.

Chemistry is molrs (CGsmiles + conformer). Architecture is molpy: a CGsmiles topology grown by ``mp.Assembler`` with
``mp.GrowthPlacer``. Packing is molpack: ``LatticeGrow`` at
2.0 Å then ``GenCanPack.with_restart`` at 2.0 Å. Hydrogen packing radius
defaults to 0.2 Å (``PEO_H_RADIUS=off`` restores ``tolerance/2``):
hydrogens relax away in the first picoseconds of MD, so making them
fight for space here only costs the heavy-atom packing.

::

    python python/examples/pack_peo_linear.py 8 8 0.5 42
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import molpy as mp
import molrs
import numpy as np
from molpy.conformer import Conformer
from molrs import Atomistic

import molpack

OUT = Path(__file__).resolve().parent / "out"

PEO_C_INF = 5.5
TET = 1.9106332


EO_UNIT = "[<]OCC[>]"  # -O-CH2-CH2-, ports on O (<) and C (>)
CORE_UNIT = "C(C[>])(C[>])(C[>])C[>]"  # pentaerythritol-like four-arm core


def _unit(name: str, body: str, seed: int) -> mp.Atomistic:
    """One CGsmiles unit with its ports, as a 3D molecule with hydrogens."""
    template = molrs.io.SmilesIR.from_fragment(body).to_template()
    return Conformer(seed=seed).generate(template)[0]


def _grow(topology: str, library: dict[str, mp.Atomistic]) -> Atomistic:
    """Grow the CGsmiles ``topology`` from ``library`` into one molecule."""
    sites = mp.CGSmilesIR(topology).to_coarsegrain()
    return mp.Assembler(library, mp.GrowthPlacer()).assemble(sites, mp.Atomistic)


def linear_topology(n: int) -> str:
    return f"{{[#EO]|{n}}}"


def make_linear(n: int, *, seed: int = 42) -> Atomistic:
    if n < 1:
        raise ValueError(f"linear n must be >= 1, got {n}")
    return _grow(linear_topology(n), {"EO": _unit("EO", EO_UNIT, seed)})


def _target(polymer: Atomistic, n_mol: int, name: str) -> molpack.Target:
    frame = polymer.to_frame()
    target = molpack.Target(frame, n_mol).with_name(name)
    raw = os.environ.get("PEO_H_RADIUS", "0.2")
    if raw not in ("", "off", "none"):
        elems = list(frame["atoms"].view("element"))
        h_idx = [i for i, e in enumerate(elems) if str(e).strip() == "H"]
        if h_idx:
            r = float(raw)
            target = target.with_atom_radius(h_idx, r)
            print(f"  H radius     : {r} Å  ({len(h_idx)} hydrogens)")
    return target


def pack_linear(n: int, n_mol: int, density: float, seed: int):
    polymer = make_linear(n, seed=seed)
    print(f"  topology     : {linear_topology(n)}")
    n_at = polymer.n_atoms
    n_bd = len(list(polymer.bonds))
    shape = "tree" if n_bd == n_at - 1 else ("unicyclic" if n_bd == n_at else "cyclic+")
    print(f"  template     : {n_at} atoms, {n_bd} bonds  ({shape})")
    print(f"  copies       : {n_mol}   density {density} g/cm³")
    target = _target(polymer, n_mol, "lin-PEO")
    progress = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "0") != "0"
    prior = molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)
    print(
        "  lattice      : LatticeGrow occupancy-guard @ 2.0 Å → "
        "GenCanPack.with_restart @ 2.0 Å"
    )
    grown = (
        molpack.LatticeGrow(prior)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_density(density)
        .with_progress(progress)
        .run([target], max_loops=max(40, n_mol * 8))
    )
    print(
        f"  grow         : converged={grown.converged}  "
        f"fdist={grown.fdist:.4e}  degraded={grown.degraded}  "
        f"intra scored {grown.intra.scored:.3f} Å"
    )
    pushed = (
        molpack.GenCanPack()
        .with_restart(grown)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_progress(progress)
        .run([target], max_loops=80)
    )
    print(f"  push-off     : converged={pushed.converged}  fdist={pushed.fdist:.4e}")
    return pushed


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    n = int(args[0]) if args else 8
    n_mol = int(args[1]) if len(args) > 1 else 8
    density = float(args[2]) if len(args) > 2 else 0.5
    seed = int(args[3]) if len(args) > 3 else 42
    print("── linear PEO (open space, molrs chemistry, molpy architecture) ──")
    t0 = time.perf_counter()
    packed = pack_linear(n, n_mol, density, seed).frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.write_mrec(str(OUT / "pack_peo_linear.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_peo_linear.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(
            str(OUT / "pack_peo_linear.dump.local"), [packed]
        )
    print(f"  wall         : {time.perf_counter() - t0:.3f} s")


if __name__ == "__main__":
    main()
