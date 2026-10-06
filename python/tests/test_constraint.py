"""Stacking molrs regions on a target — construction and attach only."""

from __future__ import annotations

import molrs
import numpy as np

import molpack


def _make_frame():
    return molrs.Frame(
        {
            "atoms": {
                "x": np.array([0.0]),
                "y": np.array([0.0]),
                "z": np.array([0.0]),
                "element": ["O"],
            }
        }
    )


class TestRegionConstruction:
    def test_box_and_sphere_reprs(self):
        box = molrs.Cuboid([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        assert "Cuboid" in repr(box)
        assert "Sphere" in repr(molrs.Sphere([0.0, 0.0, 0.0], 5.0))
        assert "composed" in repr(~molrs.Sphere([0.0, 0.0, 0.0], 2.0))

    def test_half_spaces_replace_above_and_below(self):
        below = molrs.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 10.0])
        above = ~molrs.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
        pts = np.array([[0.0, 0.0, 5.0], [0.0, 0.0, 12.0], [0.0, 0.0, -1.0]])
        assert list(below.contains(pts)) == [True, False, True]
        assert list(above.contains(pts)) == [True, True, False]


class TestConstraintStackingOnTarget:
    """Stacking multiple restraints via Target.with_restraint() calls."""

    def test_two_restraints(self):
        t = (
            molpack.Target(_make_frame(), count=1)
            .with_restraint(molrs.Cuboid([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]))
            .with_restraint(molrs.Sphere([5.0, 5.0, 5.0], 5.0))
        )
        assert t is not None

    def test_three_restraints(self):
        t = (
            molpack.Target(_make_frame(), count=1)
            .with_restraint(molrs.Cuboid([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]))
            .with_restraint(molrs.Sphere([5.0, 5.0, 5.0], 5.0))
            .with_restraint(~molrs.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]))
        )
        assert t is not None

    def test_all_region_shapes(self):
        t = molpack.Target(_make_frame(), count=1)
        for r in [
            molrs.Cuboid([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
            molrs.Sphere([0.0, 0.0, 0.0], 5.0),
            ~molrs.Sphere([0.0, 0.0, 0.0], 2.0),
            ~molrs.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
            molrs.HalfSpace([0.0, 0.0, 1.0], [0.0, 0.0, 10.0]),
            molrs.Cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 3.0, 10.0),
            molrs.Ellipsoid([0.0, 0.0, 0.0], [3.0, 4.0, 5.0]),
        ]:
            t = t.with_restraint(r)
        assert t is not None
