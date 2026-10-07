"""molpack — Packmol-grade molecular packing with Python bindings."""

from importlib.metadata import version as _dist_version

from ._protocols import Callback, Restraint
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
    GencanPack,
    IntraResidual,
    InvalidPbcBoxError,
    LatticeGrow,
    MaxIterationsError,
    NoTargetsError,
    PackError,
    Pipeline,
    ScriptJob,
    SelfSeparation,
    StageProgress,
    State,
    StepContext,
    StepReport,
    TabulatedPlane,
    TabulatedPoint,
    Target,
    TorsionPrior,
    init_thread_pool,
    load_script,
    num_threads,
    rayon_enabled,
)

# The installed wheel's version. The molrs compatibility check is the
# extension's import-time ABI handshake (`molrs_capsule::check_abi`), which runs
# when `.molpack` is imported above.
version: str = _dist_version("molcrafts-molpack")

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
    "GencanPack",
    "CbmcGrow",
    "LatticeGrow",
    "Pipeline",
    "State",
    "IntraResidual",
    "StepReport",
    "StageProgress",
    "StepContext",
    # Script loader (`.inp` input)
    "ScriptJob",
    "load_script",
    # Parallel evaluation (rayon)
    "rayon_enabled",
    "num_threads",
    "init_thread_pool",
    # Duck-type protocols
    "Callback",
    "Restraint",
    # Errors
    "PackError",
    "ConstraintsFailedError",
    "MaxIterationsError",
    "NoTargetsError",
    "EmptyMoleculeError",
    "InvalidPbcBoxError",
    "version",
]
