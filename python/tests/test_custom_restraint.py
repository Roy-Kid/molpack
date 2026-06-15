"""Tests for user-defined Python restraints (duck-typed `f`/`fg`).

A duck-typed object passed to ``Target.with_restraint`` is routed to the
**group** extension point: its ``f``/``fg`` receive ``coords`` — every copy's
``(x, y, z)`` at once — and ``fg`` returns ``(energy, [(gx, gy, gz), ...])`` with
one gradient triple per copy. (Built-in geometric restraints still take the
per-atom path; only Python objects route to the group form.)
"""

from __future__ import annotations

import molrs
import numpy as np
import pytest

import molpack


def _single_atom_frame() -> molrs.Frame:
    return molrs.Frame.from_dict(
        {
            "blocks": {
                "atoms": {
                    "x": np.array([0.0]),
                    "y": np.array([0.0]),
                    "z": np.array([0.0]),
                    "element": ["O"],
                }
            }
        }
    )


BOX_LO = [0.0, 0.0, 0.0]
BOX_HI = [40.0, 40.0, 40.0]


def _packer() -> molpack.Molpack:
    return molpack.Molpack().with_progress(False)


class InsideSpherePy:
    """Inside-sphere penalty expressed over the whole group (group contract).

    The penalty is the independent per-copy inside-sphere overshoot summed over
    every copy, so it is a decoupled (diagonal) special case of the group form:
    zero contribution when a copy is inside, else ``strength * d²`` with gradient
    ``2 * strength * (x - c) / |x - c|`` scattered onto that copy. ``coords`` is
    the list of every copy's ``(x, y, z)``; ``fg`` returns one gradient triple
    per copy. ``strength`` plays the role the built-in distribution restraints'
    ``strength``/``lambda`` multiplier does — the fixed restraint is not annealed
    with the radius schedule, so ``scale``/``scale2`` are unused.
    """

    def __init__(self, center: list[float], radius: float, strength: float = 1000.0):
        self.center = np.asarray(center, dtype=np.float64)
        self.radius = radius
        self.strength = float(strength)

    def _overshoot(self, x):
        rel = np.asarray(x, dtype=np.float64) - self.center
        dist = float(np.linalg.norm(rel))
        overshoot = dist - self.radius
        return rel, dist, overshoot

    def f(self, coords, scale, scale2):
        energy = 0.0
        for x in coords:
            _, _, overshoot = self._overshoot(x)
            if overshoot > 0.0:
                energy += self.strength * overshoot * overshoot
        return energy

    def fg(self, coords, scale, scale2):
        energy = 0.0
        grads: list[tuple[float, float, float]] = []
        for x in coords:
            rel, dist, overshoot = self._overshoot(x)
            if overshoot <= 0.0:
                grads.append((0.0, 0.0, 0.0))
                continue
            energy += self.strength * overshoot * overshoot
            # d/dx (d - R)² where d = |x - c| → 2 (d - R) * (x - c) / d
            factor = 2.0 * self.strength * overshoot / dist
            grads.append(tuple(factor * rel))
        return energy, grads


class TestAttachment:
    def test_duck_typed_object_accepted(self):
        t = molpack.Target(_single_atom_frame(), count=1).with_restraint(
            InsideSpherePy([0.0, 0.0, 0.0], 5.0)
        )
        assert t is not None

    def test_missing_fg_rejected(self):
        class OnlyF:
            def f(self, x, s, s2):
                return 0.0

        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(_single_atom_frame(), count=1).with_restraint(OnlyF())

    def test_missing_f_rejected(self):
        class OnlyFg:
            def fg(self, x, s, s2):
                return 0.0, (0.0, 0.0, 0.0)

        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(_single_atom_frame(), count=1).with_restraint(OnlyFg())


class TestPackingBehavior:
    def test_confines_atoms_inside_sphere(self):
        """A duck-typed group restraint should keep every copy within the sphere.

        The duck-typed object routes to the group path, which only runs in the
        main objective phase — so we pack a non-trivial number of copies inside
        a bounding box (a built-in per-atom restraint) and let the custom
        inside-sphere group penalty pull stragglers back in.
        """
        sphere_center = np.array([20.0, 20.0, 20.0])
        radius = 12.0
        target = (
            molpack.Target(_single_atom_frame(), count=40)
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(
                InsideSpherePy(sphere_center.tolist(), radius, strength=2000.0)
            )
        )

        result = (
            _packer()
            .with_seed(1)
            .with_tolerance(2.0)
            .pack_with_report([target], max_loops=80)
        )

        positions = np.asarray(result.positions)
        distances = np.linalg.norm(positions - sphere_center, axis=1)
        # A tolerance > 0 is needed: the soft restraint allows small
        # violations, and `frest` shows the max. Slack of a few
        # tolerance radii is generous; we're validating concept, not
        # numerical parity.
        assert np.all(distances <= radius + 4.0), (
            f"atoms escaped: max dist {distances.max():.3f} > {radius + 4.0}"
        )


class TestCallContract:
    def test_fg_is_invoked_with_whole_group_and_scales(self):
        # A non-trivial copy count forces the main objective phase, where the
        # group term runs (it is gated off during independent placement).
        count = 60
        seen: list[tuple] = []

        class Recorder:
            def f(self, coords, scale, scale2):
                return 0.0

            def fg(self, coords, scale, scale2):
                seen.append((list(coords), float(scale), float(scale2)))
                return 0.0, [(0.0, 0.0, 0.0)] * len(coords)

        target = (
            molpack.Target(_single_atom_frame(), count=count)
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(Recorder())
        )
        _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
            [target], max_loops=20
        )

        assert seen, "fg was never called"
        coords, scale, scale2 = seen[0]
        assert len(coords) == count, "fg must receive every copy at once"
        assert len(coords[0]) == 3
        assert all(isinstance(v, float) for v in coords[0])
        assert isinstance(scale, float)
        assert isinstance(scale2, float)


class TestErrorPropagation:
    def test_exception_in_fg_is_reraised(self):
        class Explodes:
            def f(self, coords, s, s2):
                return 0.0

            def fg(self, coords, s, s2):
                raise ValueError("boom from restraint")

        target = (
            molpack.Target(_single_atom_frame(), count=60)
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(Explodes())
        )
        with pytest.raises(ValueError, match="boom from restraint"):
            _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
                [target], max_loops=20
            )

    def test_fg_wrong_return_shape_is_reraised(self):
        class WrongShape:
            def f(self, coords, s, s2):
                return 0.0

            def fg(self, coords, s, s2):
                return 0.0  # missing the (energy, grads) tuple

        target = (
            molpack.Target(_single_atom_frame(), count=60)
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(WrongShape())
        )
        with pytest.raises(TypeError, match="fg.* must return"):
            _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
                [target], max_loops=20
            )
