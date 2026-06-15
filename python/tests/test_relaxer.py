"""Tests for ``molpack.relaxer.LAMMPSRelaxer`` (post-pack system relaxation).

The end-to-end relaxation needs a LAMMPS binary; those tests skip when none is
found.  The unit tests exercise the façade itself and need no external tools.
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest

import molpack
from molpack.relaxer import LAMMPSRelaxer

_LMP = next((c for c in ("lmp", "lmp_serial", "lmp_mpi") if shutil.which(c)), None)
_NEEDS_LAMMPS = pytest.mark.skipif(_LMP is None, reason="needs a LAMMPS binary")


# --------------------------------------------------------------------------
# Unit: the façade, with no LAMMPS involved
# --------------------------------------------------------------------------


def test_exposed_via_relaxer_submodule() -> None:
    """The required surface is ``molpack.relaxer.LAMMPSRelaxer``."""
    assert molpack.relaxer.LAMMPSRelaxer is LAMMPSRelaxer


def test_construction_stores_config() -> None:
    ff = object()  # opaque; the façade only holds and forwards it
    relaxer = LAMMPSRelaxer(ff, executable="lmp_mpi", launcher=["mpirun", "-np", "4"])
    assert relaxer.ff is ff
    assert relaxer.executable == "lmp_mpi"
    assert relaxer.launcher == ["mpirun", "-np", "4"]
    assert "lmp_mpi" in repr(relaxer)


def test_accepts_pack_result_or_frame() -> None:
    """A ``PackResult`` is unwrapped to its ``.frame``; a frame passes through."""
    from molpack.relaxer import _frame_of

    sentinel = object()

    class _FakeResult:
        frame = sentinel

    assert _frame_of(_FakeResult()) is sentinel
    bare = object()
    assert _frame_of(bare) is bare


# --------------------------------------------------------------------------
# End-to-end: pack with molpack, relax with LAMMPS through the façade
# --------------------------------------------------------------------------


def _packed_dimers():
    """Pack 20 neutral C-C dimers in a periodic box; return (PackResult, ff)."""
    import molrs
    from molpy.core.forcefield import AtomStyle, BondStyle, ForceField, PairStyle

    template = molrs.Frame.from_dict(
        {
            "blocks": {
                "atoms": {
                    "x": np.array([0.0, 1.5]),
                    "y": np.zeros(2),
                    "z": np.zeros(2),
                    "element": ["C", "C"],
                    "type": ["C", "C"],
                    "charge": np.zeros(2),
                },
                "bonds": {
                    "atomi": np.array([0], dtype=np.int64),
                    "atomj": np.array([1], dtype=np.int64),
                    "type": ["C-C"],
                },
            }
        }
    )
    lo, hi = [0.0, 0.0, 0.0], [22.0, 22.0, 22.0]
    target = molpack.Target(template, 20).with_restraint(
        molpack.InsideBoxRestraint(lo, hi)
    )
    result = (
        molpack.Molpack()
        .with_tolerance(2.0)
        .with_seed(1)
        .with_periodic_box(lo, hi)
        .pack_with_report([target], max_loops=50)
    )

    # `style.def_type(...)` works at runtime (def_style returns the chainable
    # ForceField); the installed molpy stub annotates def_style -> Style, which lacks
    # def_type, so ty flags it. Suppress on the def_type line until the molpy stub is
    # fixed upstream.
    ff = ForceField("cc")
    atom_style = ff.def_style(AtomStyle(name="full"))
    carbon = atom_style.def_type(  # ty: ignore[unresolved-attribute]
        "C", mass=12.011
    )
    bond_style = ff.def_style(BondStyle(name="harmonic"))
    bond_style.def_type(  # ty: ignore[unresolved-attribute]
        carbon, carbon, k=300.0, r0=1.5
    )
    pair_style = ff.def_style(PairStyle(name="lj/cut/coul/cut"))
    pair_style.def_type(  # ty: ignore[unresolved-attribute]
        carbon, carbon, epsilon=0.05, sigma=3.4
    )
    return result, ff


@_NEEDS_LAMMPS
def test_relax_pack_result_returns_relaxed_frame() -> None:
    import molrs

    result, ff = _packed_dimers()
    relaxed = LAMMPSRelaxer(ff).relax(result, max_iter=200)

    assert isinstance(relaxed, molrs.Frame)
    assert len(relaxed["atoms"].view("x")) == result.natoms  # nothing lost
    # `.box` exists at runtime; the installed molrs Frame stub omits it.
    assert relaxed.box is not None  # box preserved  # ty: ignore[unresolved-attribute]
    xyz = np.stack(
        [
            np.asarray(relaxed["atoms"].view("x")),
            np.asarray(relaxed["atoms"].view("y")),
            np.asarray(relaxed["atoms"].view("z")),
        ],
        axis=1,
    )
    assert np.isfinite(xyz).all()


@_NEEDS_LAMMPS
def test_md_settle_runs() -> None:
    result, ff = _packed_dimers()
    relaxed = LAMMPSRelaxer(ff).md(result, ensemble="nve/limit", steps=50)
    assert len(relaxed["atoms"].view("x")) == result.natoms
