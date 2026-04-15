import molpack


class TestConstraintConstruction:
    def test_inside_box(self):
        c = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        assert repr(c) == "InsideBox(...)"

    def test_inside_sphere(self):
        c = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        assert repr(c) == "InsideSphere(...)"

    def test_outside_sphere(self):
        c = molpack.OutsideSphere(2.0, [0.0, 0.0, 0.0])
        assert repr(c) == "OutsideSphere(...)"

    def test_above_plane(self):
        c = molpack.AbovePlane([0.0, 0.0, 1.0], 0.0)
        assert repr(c) == "AbovePlane(...)"

    def test_below_plane(self):
        c = molpack.BelowPlane([0.0, 0.0, 1.0], 10.0)
        assert repr(c) == "BelowPlane(...)"


class TestConstraintStackingOnTarget:
    """Stacking multiple restraints via Target.with_restraint() calls."""

    def _make_frame(self):
        import numpy as np

        return {
            "atoms": {
                "x": np.array([0.0]),
                "y": np.array([0.0]),
                "z": np.array([0.0]),
                "element": ["O"],
            }
        }

    def test_two_restraints(self):
        frame = self._make_frame()
        t = (
            molpack.Target("mol", frame, count=1)
            .with_restraint(molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]))
            .with_restraint(molpack.InsideSphere(5.0, [5.0, 5.0, 5.0]))
        )
        assert t is not None

    def test_three_restraints(self):
        frame = self._make_frame()
        t = (
            molpack.Target("mol", frame, count=1)
            .with_restraint(molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]))
            .with_restraint(molpack.InsideSphere(5.0, [5.0, 5.0, 5.0]))
            .with_restraint(molpack.AbovePlane([0.0, 0.0, 1.0], 0.0))
        )
        assert t is not None

    def test_all_restraint_types(self):
        frame = self._make_frame()
        t = molpack.Target("mol", frame, count=1)
        for r in [
            molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
            molpack.InsideSphere(5.0, [0.0, 0.0, 0.0]),
            molpack.OutsideSphere(2.0, [0.0, 0.0, 0.0]),
            molpack.AbovePlane([0.0, 0.0, 1.0], 0.0),
            molpack.BelowPlane([0.0, 0.0, 1.0], 10.0),
        ]:
            t = t.with_restraint(r)
        assert t is not None
