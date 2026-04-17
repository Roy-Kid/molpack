"""Tests for user-defined Python restraints (duck-typed `f`/`fg`)."""

from __future__ import annotations

import numpy as np
import pytest

import molpack


def _single_atom_frame() -> dict:
    return {
        "atoms": {
            "x": np.array([0.0]),
            "y": np.array([0.0]),
            "z": np.array([0.0]),
            "element": ["O"],
        }
    }


def _packer() -> molpack.Molpack:
    return molpack.Molpack(progress=False, maxit=5)


class InsideSpherePy:
    """Reimplementation of `molpack.InsideSphere` in Python.

    Penalty is quadratic in the radial overshoot, so the two-scale
    contract says we consume `scale2`. Matches
    `InsideSphereRestraint::fg` in `src/restraint.rs` branch-for-branch:
    zero contribution when inside, else `scale2 * d²` with gradient
    `2 * scale2 * (x - c)`.
    """

    def __init__(self, center: list[float], radius: float):
        self.center = np.asarray(center, dtype=np.float64)
        self.radius = radius

    def _overshoot(self, x: tuple[float, float, float]):
        rel = np.asarray(x, dtype=np.float64) - self.center
        dist = float(np.linalg.norm(rel))
        overshoot = dist - self.radius
        return rel, dist, overshoot

    def f(self, x, scale, scale2):
        _, _, overshoot = self._overshoot(x)
        if overshoot <= 0.0:
            return 0.0
        return scale2 * overshoot * overshoot

    def fg(self, x, scale, scale2):
        rel, dist, overshoot = self._overshoot(x)
        if overshoot <= 0.0:
            return 0.0, (0.0, 0.0, 0.0)
        # d/dx (d - R)² where d = |x - c| → 2 (d - R) * (x - c) / d
        factor = 2.0 * scale2 * overshoot / dist
        return scale2 * overshoot * overshoot, tuple(factor * rel)


class TestAttachment:
    def test_duck_typed_object_accepted(self):
        t = molpack.Target("mol", _single_atom_frame(), count=1).with_restraint(
            InsideSpherePy([0.0, 0.0, 0.0], 5.0)
        )
        assert t is not None

    def test_missing_fg_rejected(self):
        class OnlyF:
            def f(self, x, s, s2):
                return 0.0

        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(
                "mol", _single_atom_frame(), count=1
            ).with_restraint(OnlyF())

    def test_missing_f_rejected(self):
        class OnlyFg:
            def fg(self, x, s, s2):
                return 0.0, (0.0, 0.0, 0.0)

        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(
                "mol", _single_atom_frame(), count=1
            ).with_restraint(OnlyFg())


class TestPackingBehavior:
    def test_confines_atoms_inside_sphere(self):
        """Python restraint should keep atoms within the sphere."""
        sphere_center = np.array([0.0, 0.0, 0.0])
        radius = 6.0
        target = molpack.Target(
            "mol", _single_atom_frame(), count=10
        ).with_restraint(InsideSpherePy(sphere_center.tolist(), radius))

        result = _packer().pack([target], max_loops=80, seed=1)

        positions = np.asarray(result.positions)
        distances = np.linalg.norm(positions - sphere_center, axis=1)
        # A tolerance > 0 is needed: the soft restraint allows small
        # violations, and `frest` shows the max. Slack of a few
        # tolerance radii is generous; we're validating concept, not
        # numerical parity.
        assert np.all(distances <= radius + 2.5), (
            f"atoms escaped: max dist {distances.max():.3f} > {radius + 2.5}"
        )


class TestCallContract:
    def test_fg_is_invoked_with_three_tuple_and_scales(self):
        seen: list[tuple] = []

        class Recorder:
            def f(self, x, scale, scale2):
                return 0.0

            def fg(self, x, scale, scale2):
                seen.append((tuple(x), float(scale), float(scale2)))
                return 0.0, (0.0, 0.0, 0.0)

        target = molpack.Target(
            "mol", _single_atom_frame(), count=2
        ).with_restraint(Recorder())
        _packer().pack([target], max_loops=3, seed=1)

        assert seen, "fg was never called"
        x, scale, scale2 = seen[0]
        assert len(x) == 3
        assert all(isinstance(v, float) for v in x)
        assert isinstance(scale, float)
        assert isinstance(scale2, float)


class TestErrorPropagation:
    def test_exception_in_fg_is_reraised(self):
        class Explodes:
            def f(self, x, s, s2):
                return 0.0

            def fg(self, x, s, s2):
                raise ValueError("boom from restraint")

        target = molpack.Target(
            "mol", _single_atom_frame(), count=2
        ).with_restraint(Explodes())
        with pytest.raises(ValueError, match="boom from restraint"):
            _packer().pack([target], max_loops=3, seed=1)

    def test_fg_wrong_return_shape_is_reraised(self):
        class WrongShape:
            def f(self, x, s, s2):
                return 0.0

            def fg(self, x, s, s2):
                return 0.0  # missing the gradient tuple

        target = molpack.Target(
            "mol", _single_atom_frame(), count=2
        ).with_restraint(WrongShape())
        with pytest.raises(TypeError, match="fg.* must return"):
            _packer().pack([target], max_loops=3, seed=1)
