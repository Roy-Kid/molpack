"""molpack — Packmol-grade molecular packing with Python bindings."""

from ._protocols import Handler, Restraint
from .molpack import (
    AbovePlaneRestraint,
    Angle,
    Axis,
    BelowPlaneRestraint,
    CenteringMode,
    ConflictingPeriodicBoxesError,
    ConstraintsFailedError,
    EmptyMoleculeError,
    ExponentialPlane,
    ExponentialPoint,
    GaussianPlane,
    GaussianPoint,
    InsideBoxRestraint,
    InsideSphereRestraint,
    InvalidPBCBoxError,
    MaxIterationsError,
    Molpack,
    NoTargetsError,
    OutsideSphereRestraint,
    PackError,
    PackResult,
    ScriptJob,
    StepInfo,
    TabulatedPlane,
    TabulatedPoint,
    Target,
    init_thread_pool,
    load_script,
    num_threads,
    rayon_enabled,
)

__all__ = [
    # Typed values
    "Angle",
    "Axis",
    "CenteringMode",
    # Restraints
    "InsideBoxRestraint",
    "InsideSphereRestraint",
    "OutsideSphereRestraint",
    "AbovePlaneRestraint",
    "BelowPlaneRestraint",
    # Group-level distribution-matching restraints
    "GaussianPlane",
    "GaussianPoint",
    "ExponentialPlane",
    "ExponentialPoint",
    "TabulatedPlane",
    "TabulatedPoint",
    # Core
    "Target",
    "Molpack",
    "PackResult",
    "StepInfo",
    # Script loader (`.inp` input)
    "ScriptJob",
    "load_script",
    # Parallel evaluation (rayon)
    "rayon_enabled",
    "num_threads",
    "init_thread_pool",
    # Duck-type protocols
    "Handler",
    "Restraint",
    # Post-pack whole-system relaxation (LAMMPS via molpy)
    # Errors
    "PackError",
    "ConstraintsFailedError",
    "MaxIterationsError",
    "NoTargetsError",
    "EmptyMoleculeError",
    "InvalidPBCBoxError",
    "ConflictingPeriodicBoxesError",
]
