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
    InsideBoxRestraint,
    InsideSphereRestraint,
    InvalidPBCBoxError,
    MaxIterationsError,
    Molpack,
    NoTargetsError,
    OutsideSphereRestraint,
    PackError,
    PackResult,
    StepInfo,
    Target,
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
    # Core
    "Target",
    "Molpack",
    "PackResult",
    "StepInfo",
    # Duck-type protocols
    "Handler",
    "Restraint",
    # Errors
    "PackError",
    "ConstraintsFailedError",
    "MaxIterationsError",
    "NoTargetsError",
    "EmptyMoleculeError",
    "InvalidPBCBoxError",
    "ConflictingPeriodicBoxesError",
]
