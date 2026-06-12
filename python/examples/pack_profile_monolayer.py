"""Profile-restraint demo: a surfactant monolayer with a *measured* head profile.

Research scenario
-----------------
You have a reflectometry (neutron / X-ray) measurement of where the polar head
group of a surfactant sits relative to the water surface: a Gaussian band
centred at depth ``mu`` with width ``sigma``. You want an MD *starting*
configuration whose head groups already follow that profile, instead of the
uniform slab a hard ``inside box`` restraint would give — so the run does not
waste a long equilibration relaxing the heads into place.

Packmol's grammar can only fill hard regions uniformly; it cannot express "put
the heads on *this* curve". ``molpack``'s profile-distribution restraints can:
they bias a selected site toward a target distribution ``rho*(xi)`` by Boltzmann
inversion ``U(xi) = -kT ln(rho*/rho0)``, composing a geometric coordinate
(plane / sphere / cylinder) with a distribution (Gaussian / erf / tanh /
exponential / tabulated).

What this file shows
--------------------
The *prototype* half of the paper's extensibility story. A practitioner who
wants "heads on a Gaussian" writes the restraint in ~10 lines of Python
(:class:`PlanarGaussian` below), drops it into the *same slot* a built-in
restraint occupies (``Target.with_atom_restraint``), packs, and checks the
realised profile — no fork, no recompile. Once the idea proves out, the
identical physics is already a first-class *native* restraint, reachable from
the same Python program through the ``profile gaussian plane`` script keyword
(see :func:`native_inp` and ``profile_speed_py_vs_rust.py`` for the speed
comparison).

The Gaussian is the simplest member of the family: Boltzmann inversion of
``rho*(z) ~ exp(-(z-mu)^2 / 2 sigma^2)`` is the harmonic well
``U(z) = (kT / 2 sigma^2) (z - mu)^2`` — a spring of stiffness ``kT/sigma^2``
along the plane normal. This exactly matches molpack's native
``profile gaussian plane`` lowering (``PROFILE_KT = 1.0``), so the prototype and
the production restraint are the *same physics*.

Run
---
    cd python && maturin develop --release   # once
    python examples/pack_profile_monolayer.py
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import molrs
import numpy as np

import molpack

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "examples" / "pack_bilayer"
PALMITOYL_PDB = DATA / "palmitoil.pdb"

# palmitoil.pdb is a 42-atom palmitoyl chain. Following the canonical Packmol
# `bilayer.inp`, atoms 31 and 32 (1-based) are the polar head; here they are the
# representative site we bias onto the measured Gaussian band. Per §6.1 of the
# spec we bias a *site subset*, never every atom — biasing the whole molecule
# would distort it.
ATOMS_PER_MOL = 42
HEAD_ATOMS_0BASED = (30, 31)  # = .inp `atoms 31 32`
HEAD_ATOMS_1BASED = (31, 32)

# molpack's native `profile gaussian plane` lowering uses kT = 1.0 (an
# energy-scale factor, not a temperature). The prototype mirrors it so the two
# restraints are numerically identical.
PROFILE_KT = 1.0


class PlanarGaussian:
    """A planar Gaussian profile restraint, written as a duck-typed restraint.

    Biases a site toward ``rho*(xi) ~ exp(-(xi-mu)^2 / 2 sigma^2)`` along a
    plane normal, where ``xi = n_hat . (x - x0)`` is the signed distance to the
    reference plane. Boltzmann inversion gives the harmonic well
    ``U(xi) = (kT / 2 sigma^2) (xi - mu)^2`` with force
    ``dU/dxi = (kT / sigma^2) (xi - mu)`` along ``n_hat``.

    Implements molpack's restraint protocol — ``f`` (energy) and ``fg``
    (energy + gradient tuple) — so the packer consumes it exactly as it does a
    native restraint. The profile penalty is a *linear*-energy term, so it rides
    the global ``scale`` soft-start knob (not ``scale2``), matching the native
    ``ProfileRestraint``.
    """

    def __init__(
        self,
        normal: list[float],
        point: list[float],
        mu: float,
        sigma: float,
        kt: float = PROFILE_KT,
    ) -> None:
        n = np.asarray(normal, dtype=np.float64)
        self.n = n / np.linalg.norm(n)
        self.x0 = np.asarray(point, dtype=np.float64)
        self.mu = float(mu)
        self.k = kt / (sigma * sigma)  # spring constant kT/sigma^2

    def _xi(self, x: tuple[float, float, float]) -> float:
        return float(self.n @ (np.asarray(x, dtype=np.float64) - self.x0))

    def f(self, x, scale, scale2):  # noqa: ARG002 - scale2 unused (linear term)
        d = self._xi(x) - self.mu
        return scale * 0.5 * self.k * d * d

    def fg(self, x, scale, scale2):  # noqa: ARG002 - scale2 unused (linear term)
        d = self._xi(x) - self.mu
        energy = scale * 0.5 * self.k * d * d
        grad = scale * self.k * d * self.n  # (dU/dxi) * n_hat
        return energy, (float(grad[0]), float(grad[1]), float(grad[2]))


@dataclass(frozen=True)
class MonolayerSpec:
    """Geometry of the monolayer packing problem (shared by both paths)."""

    count: int = 64
    box_lo: tuple[float, float, float] = (0.0, 0.0, 0.0)
    box_hi: tuple[float, float, float] = (40.0, 40.0, 30.0)
    mu: float = 6.0  # measured head depth (Å)
    sigma: float = 2.0  # measured head-band width (Å)
    seed: int = 20240612
    tolerance: float = 2.0


def build_python_target(frame, spec: MonolayerSpec):
    """Build the monolayer target using the *Python* PlanarGaussian prototype."""
    head = PlanarGaussian(
        normal=[0.0, 0.0, 1.0],
        point=[0.0, 0.0, 0.0],
        mu=spec.mu,
        sigma=spec.sigma,
        kt=PROFILE_KT,
    )
    return (
        molpack.Target(frame, count=spec.count)
        .with_name("surfactant")
        .with_restraint(
            molpack.InsideBoxRestraint(list(spec.box_lo), list(spec.box_hi))
        )
        .with_atom_restraint(list(HEAD_ATOMS_0BASED), head)
    )


def native_inp(spec: MonolayerSpec, output: Path) -> str:
    """The *same* problem as a Packmol-style `.inp`, using the native keyword.

    The head bias is a single ``profile gaussian plane`` line inside the
    existing ``atoms ... end atoms`` sub-block — the production restraint that
    the Python prototype above stands in for.
    """
    lo = spec.box_lo
    hi = spec.box_hi
    a, b = HEAD_ATOMS_1BASED
    return (
        f"tolerance {spec.tolerance}\n"
        f"seed {spec.seed}\n"
        "filetype pdb\n"
        f"output {output}\n\n"
        f"structure {PALMITOYL_PDB}\n"
        f"  number {spec.count}\n"
        f"  inside box {lo[0]} {lo[1]} {lo[2]} {hi[0]} {hi[1]} {hi[2]}\n"
        f"  atoms {a} {b}\n"
        f"    profile gaussian plane 0. 0. 1. 0. 0. 0. "
        f"mu {spec.mu} sigma {spec.sigma} density\n"
        "  end atoms\n"
        "end structure\n"
    )


def head_z(positions: np.ndarray, count: int) -> np.ndarray:
    """Extract the z of every biased head site across all packed copies."""
    pos = np.asarray(positions, dtype=np.float64)
    per_mol = pos.reshape(count, ATOMS_PER_MOL, 3)
    return per_mol[:, list(HEAD_ATOMS_0BASED), 2].reshape(-1)


def summarize(zs: np.ndarray, spec: MonolayerSpec) -> str:
    within_2sigma = float(np.mean(np.abs(zs - spec.mu) <= 2.0 * spec.sigma))
    return (
        f"  target  : mu={spec.mu:.2f} Å  sigma={spec.sigma:.2f} Å\n"
        f"  realised: mean={zs.mean():.2f} Å  std={zs.std():.2f} Å  "
        f"(n={zs.size} head sites)\n"
        f"  fraction of heads within mu±2sigma: {within_2sigma:.1%}\n"
    )


def main() -> None:
    spec = MonolayerSpec()
    frame = molrs.read_pdb(str(PALMITOYL_PDB))

    show = os.environ.get("MOLPACK_EXAMPLE_PROGRESS", "1") != "0"
    target = build_python_target(frame, spec)
    packer = (
        molpack.Molpack()
        .with_progress(show)
        .with_seed(spec.seed)
        .with_tolerance(spec.tolerance)
    )
    result = packer.pack_with_report([target], max_loops=400)

    zs = head_z(result.positions, spec.count)
    print(
        f"\nMonolayer packed (Python PlanarGaussian prototype): "
        f"converged={result.converged} natoms={result.natoms} "
        f"fdist={result.fdist:.4f} frest={result.frest:.4f}"
    )
    print("Head-group depth profile vs the measured Gaussian:")
    print(summarize(zs, spec))
    print(
        "The heads track the target band (a soft restraint *approaches* rho*),\n"
        "whereas a bare `inside box` would spread them uniformly across z.\n"
        "Same physics is available natively: see native_inp() and\n"
        "profile_speed_py_vs_rust.py for the prototype-vs-production speed gap."
    )


if __name__ == "__main__":
    main()
