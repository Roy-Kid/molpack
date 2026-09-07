"""Topological PEO: 4-arm star and macrocycle, then pack.

Chemistry is molrs (SMILES + conformer, including hydrogens). Architecture
is molpy ``PolymerBuilder``. Packing is molpack. No hand-placed coordinates.

A 4-arm star needs a tetrafunctional core. Ethylene glycol only has two
reaction sites, so ``#[EO](#[EO]:n):4`` cannot branch; the graph is
``#[X4](#[EO]:n):4`` with pentaerythritol-like ``C(CO)(CO)(CO)CO``.

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

import molpy as mp
import molrs
import numpy as np
from molpy.builder.assembly import (
    MonomerLibrary,
    PolymerBuilder,
    ResiduePlacer,
    SiteMap,
    ring_cgsmiles,
    star_cgsmiles,
)
from molpy.conformer import Conformer
from molpy.core.atomistic import Atomistic

import molpack

OUT = Path(__file__).resolve().parent / "out"

ETHER = "[O;%a:1][H].[C:2][O;%b][H]>>[O:1][C:2]"
PEO_C_INF = 5.5
TET = 1.9106332
N_ARMS = 4


def smiles_3d(smiles: str, seed: int) -> Atomistic:
    """molrs SMILES parse + molrs conformer (3D, explicit H)."""
    mol, _ = Conformer(add_hydrogens=True, seed=seed).generate(
        mp.io.read_smiles(smiles)
    )
    return mol


def peo_builder(*, seed: int = 42, with_core: bool = True) -> PolymerBuilder:
    eo = smiles_3d("OCCO", seed)
    SiteMap(eo).label_elements("O", "a", "b")
    library: dict[str, Atomistic] = {"EO": eo}
    if with_core:
        core = smiles_3d("C(CO)(CO)(CO)CO", seed + 1)
        oxygens = [a for a in core.atoms if a.get("element") == "O"]
        if len(oxygens) < N_ARMS:
            raise RuntimeError("C(CO)(CO)(CO)CO must carry four hydroxyl oxygens")
        SiteMap(core).label_atoms(oxygens[:N_ARMS], *(["a"] * N_ARMS))
        library["X4"] = core
    return PolymerBuilder(
        MonomerLibrary(library),
        mp.Reaction(ETHER),
        placer=ResiduePlacer(),
    )


def make_linear(n: int, *, seed: int = 42) -> Atomistic:
    if n < 1:
        raise ValueError(f"linear n must be >= 1, got {n}")
    return peo_builder(seed=seed, with_core=False).build_linear("EO", n)


def make_star(arm_length: int, *, seed: int = 42) -> Atomistic:
    if arm_length < 1:
        raise ValueError(f"star arm_length must be >= 1, got {arm_length}")
    return peo_builder(seed=seed, with_core=True).build_star(
        "X4", "EO", n_arms=N_ARMS, arm_length=arm_length
    )


def make_ring(n: int, *, seed: int = 42) -> Atomistic:
    if n < 3:
        raise ValueError(f"ring needs n >= 3 residues, got {n}")
    return peo_builder(seed=seed, with_core=False).build_ring("EO", n)


def _h_indices(frame) -> list[int]:
    elems = list(frame["atoms"].view("element"))
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
    print(f"  topology     : {star_cgsmiles('X4', 'EO', n_arms=N_ARMS, arm_length=dp)}")
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
    print(f"  topology     : {ring_cgsmiles('EO', dp)}")
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
    print("── topological PEO (molrs chemistry, molpy architecture) ──")
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
        packed.box = molrs.Box.from_bounds(
            np.column_stack(
                [np.asarray(a["x"]), np.asarray(a["y"]), np.asarray(a["z"])]
            ),
            padding=np.ones(3),
        )
    stem = f"pack_peo_topo_{kind}"
    OUT.mkdir(parents=True, exist_ok=True)
    molrs.io.mrec.write_frame(str(OUT / f"{stem}.mrec"), packed)
    molrs.io.write_lammps_traj(str(OUT / f"{stem}.lammpstrj"), [packed])
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(str(OUT / f"{stem}.dump.local"), [packed])
    print(f"  wall         : {time.perf_counter() - t0:.3f} s")


if __name__ == "__main__":
    main()
