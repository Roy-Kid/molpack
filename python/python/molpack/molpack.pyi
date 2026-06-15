"""Type stubs for the molpack native extension (compiled from Rust via PyO3)."""

import os
from collections.abc import Sequence
from typing import Any, Protocol, Self

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Typed values
# ---------------------------------------------------------------------------

class Angle:
    ZERO: Angle
    @classmethod
    def from_degrees(cls, deg: float) -> Angle: ...
    @classmethod
    def from_radians(cls, rad: float) -> Angle: ...
    @property
    def degrees(self) -> float: ...
    @property
    def radians(self) -> float: ...
    def __eq__(self, other: object) -> bool: ...
    def __repr__(self) -> str: ...

class Axis:
    X: Axis
    Y: Axis
    Z: Axis

class CenteringMode:
    AUTO: CenteringMode
    CENTER: CenteringMode
    OFF: CenteringMode

# ---------------------------------------------------------------------------
# Restraints
# ---------------------------------------------------------------------------

class InsideBoxRestraint:
    def __init__(
        self,
        min: Sequence[float],
        max: Sequence[float],
        periodic: tuple[bool, bool, bool] = (False, False, False),
    ) -> None: ...
    def __repr__(self) -> str: ...

class InsideSphereRestraint:
    def __init__(self, center: Sequence[float], radius: float) -> None: ...
    def __repr__(self) -> str: ...

class OutsideSphereRestraint:
    def __init__(self, center: Sequence[float], radius: float) -> None: ...
    def __repr__(self) -> str: ...

class AbovePlaneRestraint:
    def __init__(self, normal: Sequence[float], distance: float) -> None: ...
    def __repr__(self) -> str: ...

class BelowPlaneRestraint:
    def __init__(self, normal: Sequence[float], distance: float) -> None: ...
    def __repr__(self) -> str: ...

# Collective (distribution-matching) restraints — drive a species' aggregate
# spatial distribution toward a target profile rather than confining each atom.
class GaussianPlane:
    def __init__(
        self,
        normal: Sequence[float],
        offset: float,
        strength: float,
        mu: float,
        sigma: float,
    ) -> None: ...
    def __repr__(self) -> str: ...

class GaussianPoint:
    def __init__(
        self, center: Sequence[float], strength: float, mu: float, sigma: float
    ) -> None: ...
    def __repr__(self) -> str: ...

class ExponentialPlane:
    def __init__(
        self,
        normal: Sequence[float],
        offset: float,
        strength: float,
        lambda_: float,
    ) -> None: ...

class ExponentialPoint:
    def __init__(
        self, center: Sequence[float], strength: float, lambda_: float
    ) -> None: ...

class TabulatedPlane:
    def __init__(
        self,
        normal: Sequence[float],
        offset: float,
        strength: float,
        xs: Sequence[float],
        rho: Sequence[float],
    ) -> None: ...

class TabulatedPoint:
    def __init__(
        self,
        center: Sequence[float],
        strength: float,
        xs: Sequence[float],
        rho: Sequence[float],
    ) -> None: ...

# Built-in restraint union — accepted wherever a native `*Restraint` is
# expected. Custom duck-typed objects (see `Restraint` Protocol in the
# package root) are also accepted.
type BuiltinRestraint = (
    InsideBoxRestraint
    | InsideSphereRestraint
    | OutsideSphereRestraint
    | AbovePlaneRestraint
    | BelowPlaneRestraint
    | GaussianPlane
    | GaussianPoint
    | ExponentialPlane
    | ExponentialPoint
    | TabulatedPlane
    | TabulatedPoint
)

class _RestraintLike(Protocol):
    def f(
        self,
        x: tuple[float, float, float],
        scale: float,
        scale2: float,
    ) -> float: ...
    def fg(
        self,
        x: tuple[float, float, float],
        scale: float,
        scale2: float,
    ) -> tuple[float, tuple[float, float, float]]: ...

type AnyRestraint = BuiltinRestraint | _RestraintLike

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------

class Target:
    def __init__(self, frame: Any, count: int) -> None: ...
    def with_name(self, name: str) -> Self: ...
    def with_restraint(self, restraint: AnyRestraint | object) -> Self: ...
    def with_atom_restraint(
        self, indices: Sequence[int], restraint: AnyRestraint | object
    ) -> Self: ...
    def with_relaxer(self, relaxer: TorsionMcRelaxer | LBFGSRelaxer) -> Self: ...
    def with_perturb_budget(self, budget: int) -> Self: ...
    def with_centering(self, mode: CenteringMode) -> Self: ...
    def with_rotation_bound(
        self, axis: Axis, center: Angle, half_width: Angle
    ) -> Self: ...
    def fixed_at(self, position: Sequence[float]) -> Self: ...
    def with_orientation(self, orientation: tuple[Angle, Angle, Angle]) -> Self: ...
    @property
    def name(self) -> str | None: ...
    @property
    def natoms(self) -> int: ...
    @property
    def count(self) -> int: ...
    @property
    def elements(self) -> list[str]: ...
    @property
    def radii(self) -> list[float]: ...
    @property
    def is_fixed(self) -> bool: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# In-loop relaxers (relaxation-assisted packing)
# ---------------------------------------------------------------------------

class TorsionMcRelaxer:
    """Monte-Carlo torsion-angle relaxer — engine-free, force-field-free.

    Built from a molecule ``frame`` (needs its bond topology); folds one chain
    during packing by proposing random rotations about its rotatable bonds and
    accepting them against the packer objective (Metropolis). Attach with
    ``Target.with_relaxer`` (requires ``count == 1``).
    """

    def __init__(self, frame: Any) -> None: ...
    def with_temperature(self, t: float) -> Self: ...
    def with_steps(self, n: int) -> Self: ...
    def with_max_delta(self, rad: float) -> Self: ...
    def with_self_avoidance(self, radius: float) -> Self: ...
    def __repr__(self) -> str: ...

class LBFGSRelaxer:
    """Force-field L-BFGS geometry relaxer (``ff`` feature).

    Built from a ``molrs``/``molpy`` force field; minimizes one molecule's
    internal geometry during packing. Attach with ``Target.with_relaxer``
    (requires ``count == 1``).
    """

    def __init__(self, forcefield: Any) -> None: ...
    def with_fmax(self, fmax: float) -> Self: ...
    def with_max_steps(self, max_steps: int) -> Self: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# StepInfo & Handler
# ---------------------------------------------------------------------------

class StepInfo:
    @property
    def loop_idx(self) -> int: ...
    @property
    def max_loops(self) -> int: ...
    @property
    def phase(self) -> int: ...
    @property
    def total_phases(self) -> int: ...
    @property
    def molecule_type(self) -> int | None: ...
    @property
    def fdist(self) -> float: ...
    @property
    def frest(self) -> float: ...
    @property
    def improvement_pct(self) -> float: ...
    @property
    def radscale(self) -> float: ...
    @property
    def precision(self) -> float: ...
    @property
    def relaxer_acceptance(self) -> list[tuple[int, float]]: ...
    def __repr__(self) -> str: ...

class _HandlerLike(Protocol):
    def on_start(self, ntotat: int, ntotmol: int) -> None: ...
    def on_step(self, info: StepInfo) -> bool | None: ...
    def on_finish(self) -> None: ...

# ---------------------------------------------------------------------------
# Molpack & PackResult
# ---------------------------------------------------------------------------

class PackResult:
    @property
    def positions(self) -> NDArray[np.float64]: ...
    @property
    def frame(self) -> Any:
        """Topology-complete ``molrs.Frame``, ready to use directly.

        Replays each target's source-frame topology (bonds/angles/dihedrals/
        impropers, with indices offset per copy) onto the packed coordinates,
        regenerates ``id`` / ``mol_id``, and stamps ``frame.box`` from the
        periodic box (if one was declared via ``with_periodic_box``). Returns a
        genuine ``molrs.Frame`` (built via the user's installed ``molrs``);
        adopt it into molpy with ``molpy.Frame.from_dict(result.frame)``. Falls
        back to a coordinates-only ``atoms`` block for ``.inp`` script packing.
        Force fields are out of scope — merge them separately.
        """
        ...
    @property
    def elements(self) -> list[str]: ...
    @property
    def converged(self) -> bool: ...
    @property
    def fdist(self) -> float: ...
    @property
    def frest(self) -> float: ...
    @property
    def natoms(self) -> int: ...
    def __repr__(self) -> str: ...

class Molpack:
    def __init__(self) -> None: ...
    def with_tolerance(self, tolerance: float) -> Self: ...
    def with_precision(self, precision: float) -> Self: ...
    def with_inner_iterations(self, n: int) -> Self: ...
    def with_init_passes(self, n: int) -> Self: ...
    def with_init_box_half_size(self, half_size: float) -> Self: ...
    def with_periodic_box(self, min: Sequence[float], max: Sequence[float]) -> Self: ...
    def with_perturb_fraction(self, f: float) -> Self: ...
    def with_random_perturb(self, enabled: bool) -> Self: ...
    def with_perturb(self, enabled: bool) -> Self: ...
    def with_avoid_overlap(self, enabled: bool) -> Self: ...
    def with_seed(self, seed: int) -> Self: ...
    def with_parallel_eval(self, enabled: bool) -> Self:
        """Run pair-kernel reductions on rayon worker threads.

        Raises :class:`RuntimeError` if ``enabled`` is true but this wheel was
        built without the ``rayon`` feature — see :func:`rayon_enabled`.
        Parallelism applies only to global evaluations; per-molecule GENCAN
        move evaluations stay serial.
        """
        ...
    def with_progress(self, enabled: bool) -> Self: ...
    def with_lammps_output(self, enabled: bool) -> Self: ...
    def with_log_level(self, level: str) -> Self: ...
    def with_log_frequency(self, n: int) -> Self: ...
    def with_handler(self, handler: _HandlerLike | object) -> Self: ...
    def with_xyz_output(self, path: str, every: int = 1) -> Self:
        """Record the packing trajectory to a multi-frame XYZ file via molpack's
        built-in ``XYZHandler`` (a frame every ``every`` loops; loop 0 included)."""
        ...
    def with_global_restraint(self, restraint: AnyRestraint | object) -> Self: ...
    def pack(self, targets: list[Target], max_loops: int = 200) -> Any: ...
    def pack_with_report(
        self, targets: list[Target], max_loops: int = 200
    ) -> PackResult: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# Script loader (`.inp` format)
# ---------------------------------------------------------------------------

class ScriptJob:
    """Output of :func:`load_script` — a packer pre-configured from a
    ``.inp`` script, plus its target list and output path.

    Supports both attribute access (``job.packer``) and tuple
    unpacking (``packer, targets, output, nloop = job``).
    """

    @property
    def packer(self) -> Molpack: ...
    @property
    def targets(self) -> list[Target]: ...
    @property
    def output(self) -> str: ...
    @property
    def nloop(self) -> int: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx: int) -> Any: ...
    def __repr__(self) -> str: ...

def load_script(path: str | os.PathLike[str]) -> ScriptJob:
    """Parse a molpack ``.inp`` script and lower it to a ready-to-run
    packer and target list.

    Relative paths inside the script (structure files, output) are
    resolved against the script's parent directory.
    """

# ---------------------------------------------------------------------------
# Parallel evaluation (rayon)
# ---------------------------------------------------------------------------

def rayon_enabled() -> bool:
    """True if this wheel was built with the ``rayon`` feature (parallel
    evaluation compiled in)."""

def num_threads() -> int:
    """Worker threads the parallel evaluator will use — the rayon global pool
    size, or ``1`` for a serial build."""

def init_thread_pool(n: int) -> None:
    """Pin the rayon global thread pool to ``n`` workers.

    Must be called before the first parallel pack — the pool is immutable once
    built. Raises :class:`RuntimeError` without the ``rayon`` feature or if the
    pool was already initialized, and :class:`ValueError` if ``n < 1``. For a
    scaling sweep, set the count once per process and launch one process per
    data point.
    """

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class PackError(RuntimeError): ...
class ConstraintsFailedError(PackError): ...
class MaxIterationsError(PackError): ...
class NoTargetsError(PackError): ...
class EmptyMoleculeError(PackError): ...
class InvalidPBCBoxError(PackError): ...
class ConflictingPeriodicBoxesError(PackError): ...
