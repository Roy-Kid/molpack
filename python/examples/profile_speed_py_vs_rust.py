"""Prototype-vs-production speed: the same profile restraint, Python vs Rust.

This is the measurement half of the extensibility story in
``pack_profile_monolayer.py``. It reuses the *identical* duck-typed
:class:`PlanarGaussian` restraint, run two ways on the same problem, both driven
from the same Python program:

  (A) Python prototype   the :class:`PlanarGaussian` restraint. The objective
                         calls back into Python once per biased site on every
                         energy / gradient evaluation (GIL + marshalling).

  (B) Native (Rust)      the same physics as the built-in ``profile gaussian
                         plane`` keyword, reached through ``molpack.load_script``
                         — the restraint runs inline in the compiled hot loop.

System: an electrode-style **diffuse layer of single-atom ions**, each biased
toward a Gaussian band at distance ``mu`` from a plane. The single-atom ion is
deliberate: when the biased site *is* the whole molecule, the restraint is
evaluated for every atom, so the per-call cost is not diluted by intramolecular
overlap work the way it is for a 2-of-42-atom head site (see the monolayer demo,
where the gap is small precisely because the restraint is a sliver of the work).
This isolates the restraint-evaluation cost and shows the true prototype-to-
production gap.

Both paths share structure, count, box, mu/sigma, seed, tolerance, and
max_loops, so the only difference is *where the restraint executes*. We report
the median wall-clock pack time of each and the Rust speed-up, plus a parity
check that the two realise the same ion profile (functional, not bitwise).

Run
---
    cd python && maturin develop --release   # once
    python examples/profile_speed_py_vs_rust.py
    # smaller/faster:
    MOLPACK_BENCH_COUNTS=200,500 MOLPACK_BENCH_LOOPS=40 MOLPACK_BENCH_RUNS=3 \
        python examples/profile_speed_py_vs_rust.py
"""

from __future__ import annotations

import os
import statistics
import tempfile
import time
from pathlib import Path

import molrs
import numpy as np

import molpack
from pack_profile_monolayer import PROFILE_KT, PlanarGaussian


def _env_ints(name: str, default: list[int]) -> list[int]:
    raw = os.environ.get(name)
    return [int(t) for t in raw.split(",")] if raw else default


COUNTS = _env_ints("MOLPACK_BENCH_COUNTS", [250, 500, 1000])
MAX_LOOPS = int(os.environ.get("MOLPACK_BENCH_LOOPS", "50"))
RUNS = int(os.environ.get("MOLPACK_BENCH_RUNS", "3"))

# Electrolyte geometry: a moderately dense slab so overlap resolution actually
# runs for the full max_loops (a trivially dilute box converges in a few steps
# and times nothing), while the per-atom restraint stays the dominant
# per-evaluation cost. Gaussian counter-ion band at mu, width sigma, along +z
# from the z=0 electrode plane.
BOX = (25.0, 25.0, 40.0)
MU = 20.0
SIGMA = 4.0
SEED = 20240612
TOLERANCE = 2.0

ION_PDB_BODY = (
    "HETATM    1 NA   ION     1       0.000   0.000   0.000"
    "  1.00  0.00          NA\n"
    "END\n"
)


def write_ion_pdb(workdir: Path) -> Path:
    path = workdir / "ion.pdb"
    path.write_text(ION_PDB_BODY)
    return path


def python_restraint() -> PlanarGaussian:
    return PlanarGaussian([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], MU, SIGMA, PROFILE_KT)


def time_python(frame, count: int) -> tuple[float, np.ndarray]:
    """Pack with the Python-callback restraint; return (seconds, ion zs)."""
    target = (
        molpack.Target(frame, count=count)
        .with_name("ion")
        .with_restraint(molpack.InsideBoxRestraint([0.0, 0.0, 0.0], list(BOX)))
        .with_restraint(python_restraint())
    )
    packer = (
        molpack.Molpack().with_progress(False).with_seed(SEED).with_tolerance(TOLERANCE)
    )
    t0 = time.perf_counter()
    result = packer.pack_with_report([target], max_loops=MAX_LOOPS)
    dt = time.perf_counter() - t0
    return dt, np.asarray(result.positions, dtype=np.float64)[:, 2]


def native_inp(ion_pdb: Path, count: int, output: Path) -> str:
    return (
        f"tolerance {TOLERANCE}\n"
        f"seed {SEED}\n"
        "filetype pdb\n"
        f"output {output}\n\n"
        f"structure {ion_pdb}\n"
        f"  number {count}\n"
        f"  inside box 0 0 0 {BOX[0]} {BOX[1]} {BOX[2]}\n"
        f"  profile gaussian plane 0. 0. 1. 0. 0. 0. mu {MU} sigma {SIGMA} density\n"
        "end structure\n"
    )


def time_native(ion_pdb: Path, count: int, workdir: Path) -> tuple[float, np.ndarray]:
    """Pack with the native `profile` keyword via load_script; return (s, zs)."""
    inp = workdir / f"ions_{count}.inp"
    out = workdir / f"ions_{count}.pdb"
    inp.write_text(native_inp(ion_pdb, count, out))

    job = molpack.load_script(inp)
    packer = job.packer.with_progress(False)
    t0 = time.perf_counter()
    result = packer.pack_with_report(job.targets, max_loops=MAX_LOOPS)
    dt = time.perf_counter() - t0
    return dt, np.asarray(result.positions, dtype=np.float64)[:, 2]


def median_time(fn, *args) -> tuple[float, np.ndarray]:
    times: list[float] = []
    zs = np.empty(0)
    for _ in range(RUNS):
        dt, zs = fn(*args)
        times.append(dt)
    return statistics.median(times), zs


def main() -> None:
    print(
        "Profile-restraint speed: Python prototype vs native Rust keyword\n"
        "  system    : single-atom ions in a diffuse layer (restraint = every atom)\n"
        f"  bias      : profile gaussian plane  mu={MU} sigma={SIGMA}  box={BOX}\n"
        f"  max_loops : {MAX_LOOPS}   runs/median : {RUNS}\n"
    )
    header = (
        f"{'ions':>6} | {'Python (s)':>11} {'Rust (s)':>9} {'speed-up':>9} "
        f"| {'parity (ion mean±std, Å)':>30}"
    )
    print(header)
    print("-" * len(header))

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        ion_pdb = write_ion_pdb(workdir)
        frame = molrs.read_pdb(str(ion_pdb))
        for count in COUNTS:
            py_t, py_z = median_time(time_python, frame, count)
            rs_t, rs_z = median_time(time_native, ion_pdb, count, workdir)
            speedup = py_t / rs_t if rs_t > 0 else float("inf")
            parity = (
                f"py {py_z.mean():.2f}±{py_z.std():.2f} / "
                f"rs {rs_z.mean():.2f}±{rs_z.std():.2f}"
            )
            print(
                f"{count:>6} | {py_t:>11.3f} {rs_t:>9.3f} {speedup:>8.1f}× "
                f"| {parity:>30}"
            )

    print(
        "\nSame restraint, same problem: the Python prototype pays a per-site\n"
        "Python call on every objective/gradient evaluation; the native keyword\n"
        "runs inline. Prototype in Python to find the physics, sink to Rust for\n"
        "the production run — the driver script is unchanged."
    )


if __name__ == "__main__":
    main()
