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
    ScriptJob,
    StepInfo,
    Target,
    load_script,
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
    # Script loader (`.inp` input)
    "ScriptJob",
    "load_script",
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
