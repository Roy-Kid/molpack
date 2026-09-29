"""molpack — Packmol-grade molecular packing with Python bindings."""

from ._protocols import Handler, Restraint
from .molpack import (
    Angle,
    AnglePrior,
    Axis,
    CbmcGrow,
    CenteringMode,
    ConstraintsFailedError,
    EmptyMoleculeError,
    ExponentialPlane,
    ExponentialPoint,
    GaussianPlane,
    GaussianPoint,
    GenCanPack,
    IntraResidual,
    InvalidPBCBoxError,
    LatticeGrow,
    MaxIterationsError,
    NoTargetsError,
    PackError,
    Pipeline,
    ScriptJob,
    SelfSeparation,
    StageInfo,
    State,
    StepContext,
    StepInfo,
    TabulatedPlane,
    TabulatedPoint,
    Target,
    TorsionPrior,
    init_thread_pool,
    load_script,
    num_threads,
    rayon_enabled,
)
from .version import MOLRS_MINOR, check_molrs_version, version

__all__ = [
    # Typed values
    "Angle",
    "Axis",
    "CenteringMode",
    # Growth statistics inputs
    "TorsionPrior",
    "AnglePrior",
    # Group-level distribution-matching restraints
    "GaussianPlane",
    "GaussianPoint",
    "ExponentialPlane",
    "ExponentialPoint",
    "TabulatedPlane",
    "TabulatedPoint",
    # Group-level separation restraint
    "SelfSeparation",
    # Core
    "Target",
    "GenCanPack",
    "CbmcGrow",
    "LatticeGrow",
    "Pipeline",
    "State",
    "IntraResidual",
    "StepInfo",
    "StageInfo",
    "StepContext",
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
    # Errors
    "PackError",
    "ConstraintsFailedError",
    "MaxIterationsError",
    "NoTargetsError",
    "EmptyMoleculeError",
    "InvalidPBCBoxError",
    "MOLRS_MINOR",
    "check_molrs_version",
    "version",
]
