"""Post-pack, whole-system relaxation via an external MD engine.

The relaxers exposed at the top level (:class:`molpack.LBFGSRelaxer`,
:class:`molpack.TorsionMcRelaxer`) are *in-loop, per-molecule* relaxers: they
reshape a single molecule's reference geometry **during** packing.  This module
is a different axis — it takes the **finished** packed box and hands the whole
assembled frame to an external engine (LAMMPS) for an energy minimisation or a
short MD settle, returning the relaxed frame.

The work is delegated to :class:`molpy.engine.LAMMPS` (``molcrafts-molpy``,
a hard dependency of molpack).

Example::

    import molpack
    from molpack.relaxer import LAMMPSRelaxer

    result = molpack.Molpack().with_periodic_box(lo, hi).pack_with_report(targets)
    relaxer = LAMMPSRelaxer(ff)               # ff: a typified molpy ForceField
    relaxed = relaxer.relax(result)           # -> molrs.Frame, overlaps minimised
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from molpy.engine import LAMMPS

if TYPE_CHECKING:
    import molrs

__all__ = ["LAMMPSRelaxer"]


def _frame_of(target: Any) -> molrs.Frame:
    """Return the frame carried by *target* (a ``PackResult`` or a ``Frame``)."""
    return getattr(target, "frame", target)


class LAMMPSRelaxer:
    """Relax a packed system with LAMMPS (minimisation or short MD).

    A thin, force-field-bound façade over :class:`molpy.engine.LAMMPS`.  The
    force field is fixed at construction; each call accepts a packed result (or
    a bare frame) and returns a new relaxed :class:`~molrs.Frame` — the input is
    never mutated.

    Args:
        ff: A typified force field (molpy ``ForceField``) supplying the
            pair/bond/angle/... coefficients.
        executable: LAMMPS binary; ``None`` auto-detects ``lmp`` / ``lmp_serial``
            / ``lmp_mpi`` on ``PATH``.
        launcher: Optional MPI / scheduler prefix, e.g. ``["mpirun", "-np", "8"]``.
        pair_style: Default LAMMPS ``pair_style`` line used for relaxation.
        atom_style: LAMMPS ``atom_style`` (``"full"`` by default).
        units: LAMMPS ``units`` (``"real"`` by default).
        workdir: Default directory for input/output files; a temporary directory
            is used per call when ``None``.
    """

    def __init__(
        self,
        ff: Any,
        *,
        executable: str | None = None,
        launcher: list[str] | None = None,
        pair_style: str = "lj/cut/coul/cut 10.0",
        atom_style: str = "full",
        units: str = "real",
        workdir: str | Path | None = None,
    ) -> None:
        self.ff = ff
        self.executable = executable
        self.launcher = launcher
        self.pair_style = pair_style
        self.atom_style = atom_style
        self.units = units
        self.workdir = workdir

    def minimize(self, target: Any, **options: Any) -> molrs.Frame:
        """Energy-minimise *target* and return the relaxed frame.

        Args:
            target: A ``PackResult`` (its ``.frame`` is used) or a ``molrs.Frame``.
                The frame must carry a periodic box.
            **options: Forwarded to :meth:`molpy.engine.LAMMPS.minimize`
                (``etol``, ``ftol``, ``max_iter``, ``max_eval``, ``pair_style``,
                ``workdir``, ...).
        """
        return self._run("minimize", target, options)

    def md(self, target: Any, **options: Any) -> molrs.Frame:
        """Run short MD on *target* (settle) and return the resulting frame.

        Args:
            target: A ``PackResult`` or a ``molrs.Frame`` with a periodic box.
            **options: Forwarded to :meth:`molpy.engine.LAMMPS.md`
                (``ensemble``, ``steps``, ``temperature``, ``timestep``, ...).
        """
        return self._run("md", target, options)

    def relax(self, target: Any, **options: Any) -> molrs.Frame:
        """Default post-pack relaxation — an alias for :meth:`minimize`."""
        return self.minimize(target, **options)

    __call__ = relax

    def _run(self, mode: str, target: Any, options: dict[str, Any]) -> molrs.Frame:
        """Resolve the engine, apply defaults, and dispatch to *mode*."""
        engine = self._engine()
        options.setdefault("pair_style", self.pair_style)
        options.setdefault("atom_style", self.atom_style)
        options.setdefault("units", self.units)
        if self.workdir is not None:
            options.setdefault("workdir", self.workdir)
        return getattr(engine, mode)(_frame_of(target), self.ff, **options)

    def _engine(self) -> Any:
        """Construct the backing :class:`molpy.engine.LAMMPS`."""
        return LAMMPS(self.executable, launcher=self.launcher)

    def __repr__(self) -> str:
        binary = self.executable or "auto"
        return f"<LAMMPSRelaxer executable='{binary}', pair_style='{self.pair_style}'>"
