"""Topological PEO: 4-arm star and macrocycle, then pack.

Chemistry and architecture are molrs: SMILES + conformer (hydrogens
included) for each unit, and CGsmiles topologies grown by
``molrs.builder.Assembler`` with ``molrs.builder.GrowthPlacer``. Packing is
molpack. No hand-placed coordinates.

A 4-arm star needs a tetrafunctional core. An EO unit only has two ports,
so the core is ``X4`` (``C(C[>])(C[>])(C[>])C[>]``, pentaerythritol-like)
with four EO arms.

::

    python python/examples/pack_peo_topo.py star 4 8 0.5 42
    python python/examples/pack_peo_topo.py ring 6 8 0.4 42

Star packing is an explicit pick: ``LatticeGrow`` at 2.0 Å then
``GenCanPack.with_restart`` at 2.0 Å. Hydrogen packing radius defaults
to 0.2 Å (``PEO_H_RADIUS=off`` restores ``tolerance/2``):
hydrogens relax away in the first picoseconds of MD, so making them
fight for space here only costs the heavy-atom packing.
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


def ring_topology(n: int) -> str:
    return "{[#EO]1" + "[#EO]" * (n - 2) + "[#EO]1}"


def make_ring(n: int, *, seed: int = 42) -> Atomistic:
    if n < 3:
        raise ValueError(f"ring needs n >= 3 residues, got {n}")
    return _grow(ring_topology(n), {"EO": _unit("EO", EO_UNIT, seed)})


def _h_indices(frame) -> list[int]:
    elems = list(frame["atoms"]["element"])
    return [i for i, e in enumerate(elems) if str(e).strip() == "H"]


def _target(polymer: Atomistic, n_mol: int, name: str) -> molpack.Target:
    frame = polymer.to_frame()
    target = molpack.Target(frame, n_mol).with_name(name)
    raw = os.environ.get("PEO_H_RADIUS", "0.2")
    if raw not in ("", "off", "none"):
        h_idx = _h_indices(frame)
        if h_idx:
            r = float(raw)
            target = target.with_atom_radius(h_idx, r)
            print(f"  H radius     : {r} Å  ({len(h_idx)} hydrogens)")
    return target


def _prior() -> molpack.TorsionPrior:
    return molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)


def _progress() -> bool:
    return os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "0") != "0"


def _report_graph(polymer: Atomistic) -> tuple[int, int]:
    n_at = polymer.n_atoms
    n_bd = len(list(polymer.bonds))
    shape = "tree" if n_bd == n_at - 1 else ("unicyclic" if n_bd == n_at else "cyclic+")
    print(f"  template     : {n_at} atoms, {n_bd} bonds  ({shape})")
    return n_at, n_bd


def lattice_then_push(
    targets: list[molpack.Target],
    *,
    density: float,
    seed: int,
    max_grow: int,
    max_push: int = 80,
) -> molpack.State:
    """Open-space melt: LatticeGrow occupancy-guard @ 2.0 Å, then GENCAN push-off."""
    prior = _prior()
    print(
        "  lattice      : LatticeGrow occupancy-guard @ 2.0 Å → GenCanPack.with_restart @ 2.0 Å"
    )
    grown = (
        molpack.LatticeGrow(prior)
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_density(density)
        .with_progress(_progress())
        .run(targets, max_loops=max_grow)
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
        .with_progress(_progress())
        .run(targets, max_loops=max_push)
    )
    print(f"  push-off     : converged={pushed.converged}  fdist={pushed.fdist:.4e}")
    return pushed


def pack_star(dp: int, n_mol: int, density: float, seed: int) -> molpack.State:
    polymer = make_star(dp, seed=seed)
    print(f"  topology     : {star_topology(dp)}")
    _report_graph(polymer)
    print(f"  copies       : {n_mol}   density {density} g/cm³")
    target = _target(polymer, n_mol, "star-PEO")
    return lattice_then_push(
        [target],
        density=density,
        seed=seed,
        max_grow=max(40, n_mol * 8),
    )


def pack_ring(dp: int, n_mol: int, density: float, seed: int) -> molpack.State:
    polymer = make_ring(dp, seed=seed)
    print(f"  topology     : {ring_topology(dp)}")
    _report_graph(polymer)
    print(f"  copies       : {n_mol}   density {density} g/cm³")
    target = _target(polymer, n_mol, "c-PEO")
    prior = _prior()
    print("  grow         : CbmcGrow (expect RingTemplate)")
    try:
        molpack.CbmcGrow(prior).with_seed(seed).with_tolerance(2.0).with_density(
            density
        ).with_progress(False).run([target], max_loops=4)
        raise RuntimeError("CbmcGrow accepted a ring — RingTemplate should have fired")
    except ValueError as err:
        if "ring" not in str(err).lower():
            raise
        print(f"  grow         : named reject — {err}")
    packed = (
        molpack.GenCanPack()
        .with_seed(seed)
        .with_tolerance(2.0)
        .with_density(density)
        .with_progress(_progress())
        .run([target], max_loops=max(40, n_mol * 8))
    )
    print(
        f"  pack         : converged={packed.converged}  "
        f"fdist={packed.fdist:.4e}  n={packed.natoms}"
    )
    return packed


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    kind = args[0] if args else "star"
    dp = int(args[1]) if len(args) > 1 else 4
    n_mol = int(args[2]) if len(args) > 2 else 8
    density = float(args[3]) if len(args) > 3 else 0.5
    seed = int(args[4]) if len(args) > 4 else 42
    print("── topological PEO (molrs chemistry and architecture) ──")
    t0 = time.perf_counter()
    if kind == "ring":
        state = pack_ring(dp, n_mol, density, seed)
    elif kind == "star":
        state = pack_star(dp, n_mol, density, seed)
    else:
        raise SystemExit("usage: pack_peo_topo.py <star|ring> [dp n_mol density seed]")
    packed = state.frame
    if packed.box is None:
        a = packed["atoms"]
        packed.box = molrs.core.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    stem = f"pack_peo_topo_{kind}"
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.write_mrec_frame(str(OUT / f"{stem}.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / f"{stem}.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(str(OUT / f"{stem}.dump.local"), [packed])
    print(f"  wall         : {time.perf_counter() - t0:.3f} s")


if __name__ == "__main__":
    main()
