"""
Type stubs for the molpack native extension (compiled from Rust via PyO3).
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Constraint types
# ---------------------------------------------------------------------------

AnyConstraint = "InsideBox | InsideSphere | OutsideSphere | AbovePlane | BelowPlane | MoleculeConstraint"

class InsideBox:
    def __init__(self, min: Sequence[float], max: Sequence[float]) -> None: ...
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

class InsideSphere:
    def __init__(self, radius: float, center: Sequence[float]) -> None: ...
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

class OutsideSphere:
    def __init__(self, radius: float, center: Sequence[float]) -> None: ...
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

class AbovePlane:
    def __init__(self, normal: Sequence[float], distance: float) -> None: ...
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

class BelowPlane:
    def __init__(self, normal: Sequence[float], distance: float) -> None: ...
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

class MoleculeConstraint:
    def and_(
        self,
        other: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> MoleculeConstraint: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------

class Target:
    @staticmethod
    def from_coords(
        positions: NDArray[np.float64],
        radii: NDArray[np.float64],
        count: int,
        elements: list[str] | None = None,
    ) -> Target: ...
    def with_name(self, name: str) -> Target: ...
    def with_constraint(
        self,
        constraint: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> Target: ...
    def with_constraint_for_atoms(
        self,
        indices: list[int],
        constraint: InsideBox
        | InsideSphere
        | OutsideSphere
        | AbovePlane
        | BelowPlane
        | MoleculeConstraint,
    ) -> Target: ...
    def with_maxmove(self, maxmove: int) -> Target: ...
    def with_center(self) -> Target: ...
    def without_centering(self) -> Target: ...
    def constrain_rotation_x(
        self, center_deg: float, half_width_deg: float
    ) -> Target: ...
    def constrain_rotation_y(
        self, center_deg: float, half_width_deg: float
    ) -> Target: ...
    def constrain_rotation_z(
        self, center_deg: float, half_width_deg: float
    ) -> Target: ...
    def fixed_at(self, position: Sequence[float]) -> Target: ...
    def fixed_at_with_euler(
        self, position: Sequence[float], euler: Sequence[float]
    ) -> Target: ...
    @property
    def natoms(self) -> int: ...
    @property
    def count(self) -> int: ...
    @property
    def elements(self) -> list[str]: ...
    @property
    def is_fixed(self) -> bool: ...
    def __repr__(self) -> str: ...

# ---------------------------------------------------------------------------
# Packer & PackResult
# ---------------------------------------------------------------------------

class PackResult:
    @property
    def positions(self) -> NDArray[np.float64]: ...
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

class Packer:
    def __init__(self, tolerance: float = 2.0, precision: float = 0.01) -> None: ...
    def with_tolerance(self, tolerance: float) -> Packer: ...
    def with_precision(self, precision: float) -> Packer: ...
    def with_maxit(self, maxit: int) -> Packer: ...
    def with_nloop0(self, nloop0: int) -> Packer: ...
    def with_sidemax(self, sidemax: float) -> Packer: ...
    def with_movefrac(self, movefrac: float) -> Packer: ...
    def with_movebadrandom(self, enabled: bool) -> Packer: ...
    def with_disable_movebad(self, disabled: bool) -> Packer: ...
    def with_pbc(self, min: Sequence[float], max: Sequence[float]) -> Packer: ...
    def with_pbc_box(self, lengths: Sequence[float]) -> Packer: ...
    def with_progress(self, enabled: bool) -> Packer: ...
    def pack(
        self,
        targets: list[Target],
        max_loops: int = 200,
        seed: int | None = None,
    ) -> PackResult: ...
    def __repr__(self) -> str: ...
