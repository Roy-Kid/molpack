"""Linear PEO grown through the solvent-accessible void of a bead-spring frame.

The frame is a LAMMPS data file: a coarse-grained polymer in a solvent, in
LJ units. Every atom that appears in ``Bonds`` is the polymer; everything
else is solvent and is dropped — the void is what the solvent occupied.

The region is built from the atoms, in memory, with no mesh: one sphere per
polymer bead of radius ``bead radius + probe radius`` (``molrs.SphereUnion``,
minimum image on the box's periodic axes) is the solvent-accessible volume,
and ``~polymer`` is the space a PEO atom centre may occupy. ``LatticeGrow``
walks Region ∩ lattice, so the chains thread the solvent channels by
construction.

Units: the file is in σ; ``SIGMA_A`` converts to Å (4.18 Å/σ by default), and
``PROBE_A`` is the PEO atom packing radius (half of ``tolerance``). The scene
stops after the grow — at melt density the residual contacts belong to the
force field downstream, and a rigid push-off inside a void has nowhere to put
them except through the beads.

::

    python python/examples/pack_peo_void.py <frame.data> 25 200 42
    python python/examples/pack_peo_void.py <frame.data> 25 200 42 4.18 1.0
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import molpy as mp
import molrs
import numpy as np
from molpy.conformer import Conformer
from molrs import Atomistic

import molpack

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
PEO_C_INF = 5.5
TET = 1.9106332
#: Å per LJ σ.
SIGMA_A = 4.18
#: Bead radius in σ (LJ contact is 2^(1/6) σ; half a σ is the hard core).
BEAD_RADIUS_SIGMA = 0.5
#: PEO atom packing radius, Å — half of the packing tolerance.
PROBE_A = 1.0


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


def select_polymer(frame) -> np.ndarray:
    """Coordinates (N, 3) of every atom that appears in a bond, file units.

    molrs's LAMMPS reader stores ``bonds["atomi"]`` / ``["atomj"]`` as
    0-based row indices into the atoms block (the file's atom ids, which may
    be permuted, are already resolved), so the rows are indexed directly.
    """
    atoms = frame["atoms"]
    bonds = frame["bonds"]
    rows = np.unique(
        np.concatenate([np.asarray(bonds["atomi"]), np.asarray(bonds["atomj"])])
    )
    xyz = np.column_stack(
        [np.asarray(atoms["x"]), np.asarray(atoms["y"]), np.asarray(atoms["z"])]
    )
    return xyz[rows]


def void_region(frame, sigma: float = SIGMA_A, probe: float = PROBE_A):
    """``~SphereUnion`` over the polymer beads, in Å, periodic like the frame's box."""
    centers = select_polymer(frame) * sigma
    h = np.asarray(frame.box.h) * sigma
    origin = np.asarray(frame.box.origin) * sigma
    box = molrs.Box(h, origin, np.asarray(frame.box.pbc))
    polymer = molrs.SphereUnion(centers, BEAD_RADIUS_SIGMA * sigma + probe, box=box)
    return ~polymer, box


def main(argv: list[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        raise SystemExit(__doc__)
    data = Path(args[0])
    n = int(args[1]) if len(args) > 1 else 25
    n_mol = int(args[2]) if len(args) > 2 else 200
    seed = int(args[3]) if len(args) > 3 else 42
    sigma = float(args[4]) if len(args) > 4 else SIGMA_A
    probe = float(args[5]) if len(args) > 5 else PROBE_A

    t0 = time.perf_counter()
    frame = molrs.io.read_lammps_data(str(data))
    void, box = void_region(frame, sigma, probe)
    lo = np.asarray(box.origin)
    hi = lo + np.diag(np.asarray(box.h))
    grid = np.stack(
        np.meshgrid(
            *[np.linspace(a, b, 30, endpoint=False) for a, b in zip(lo, hi)],
            indexing="ij",
        ),
        -1,
    ).reshape(-1, 3)
    print(
        "── PEO through the solvent-accessible void (LatticeGrow inside ~SphereUnion) ──"
    )
    print(
        f"  frame        : {data.name}  {frame['atoms'].nrows} atoms, {select_polymer(frame).shape[0]} polymer beads"
    )
    print(
        f"  scale        : {sigma:g} Å/σ   sphere radius {BEAD_RADIUS_SIGMA * sigma + probe:.2f} Å"
    )
    print(f"  cell         : [{lo.round(2).tolist()}, {hi.round(2).tolist()}] Å")
    print(f"  void fraction: {void.contains(grid).mean():.3f} on a 30³ grid")

    template = make_linear(n, seed=seed).to_frame()
    hs = [
        i for i, e in enumerate(template["atoms"]["element"]) if str(e).strip() == "H"
    ]
    target = (
        molpack.Target(template, n_mol)
        .with_name("lin-PEO")
        .with_restraint(void)
        .with_atom_radius(hs, 0.2)
    )
    prior = molpack.TorsionPrior.three_state_from_c_inf(PEO_C_INF, TET)
    grown = (
        molpack.LatticeGrow(prior)
        .with_seed(seed)
        .with_tolerance(2.0 * probe)
        .with_periodic_box(lo.tolist(), hi.tolist())
        .run([target], max_loops=max(40, n_mol * 8))
    )
    print(
        f"  grow         : converged={grown.converged}  fdist={grown.fdist:.4e}  "
        f"frest={grown.frest:.4e}  degraded={grown.degraded}"
    )
    inside = void.contains(grown.positions)
    print(f"  atoms in void: {inside.mean() * 100:.2f}%")

    OUT.mkdir(parents=True, exist_ok=True)
    packed = grown.frame
    molrs.io.write_mrec(str(OUT / "pack_peo_void.mrec"), packed)
    molrs.io.write_lammps_trajectory(
        str(OUT / "pack_peo_void.lammpstrj"),
        [packed],
        columns=["id", "element", "mol", "x", "y", "z"],
    )
    if "bonds" in packed and packed["bonds"].nrows:
        molrs.io.write_lammps_dump_local(
            str(OUT / "pack_peo_void.dump.local"), [packed]
        )
    print(f"  wall         : {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
