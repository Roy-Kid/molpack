"""Runtime-checkable duck-type protocols.

These describe the contracts that user-defined handlers and restraints
must satisfy. They are pure Python and never passed into Rust — the
native extension extracts methods by name.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from .molpack import StepInfo


@runtime_checkable
class Handler(Protocol):
    """Progress-handler contract.

    All methods are optional. Missing ones are silently skipped by the
    native extension. Returning ``True`` from :meth:`on_step` requests
    early termination.
    """

    def on_start(self, ntotat: int, ntotmol: int) -> None: ...

    def on_step(self, info: StepInfo) -> bool | None: ...

    def on_finish(self) -> None: ...


@runtime_checkable
class Restraint(Protocol):
    """Custom-restraint contract — the group-level extension point.

    A duck-typed restraint sees **every copy** of a species at once: ``coords``
    is a sequence of ``(x, y, z)`` tuples (Å), one per copy. This lets the
    penalty depend on the *joint* configuration — e.g. matching the empirical
    density of the group to a target distribution — so its gradient may couple
    the copies together. ``scale`` is the linear scaling factor (≈ distance
    tolerance); ``scale2`` is its square — both are supplied so the restraint
    can pick whichever matches the penalty's polynomial degree (linear vs
    quadratic overshoot).

    ``fg`` returns ``(energy, grads)`` where ``grads`` has the same length as
    ``coords`` and ``grads[i]`` is the ``(gx, gy, gz)`` gradient of the penalty
    with respect to ``coords[i]``.
    """

    def f(
        self,
        coords: list[tuple[float, float, float]],
        scale: float,
        scale2: float,
    ) -> float: ...

    def fg(
        self,
        coords: list[tuple[float, float, float]],
        scale: float,
        scale2: float,
    ) -> tuple[float, list[tuple[float, float, float]]]: ...
