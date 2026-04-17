import molpack


class TestConstraintConstruction:
    def test_inside_box(self):
        c = molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        assert "InsideBoxRestraint" in repr(c)
        assert "periodic=(false, false, false)" in repr(c)

    def test_inside_box_periodic(self):
        c = molpack.InsideBoxRestraint(
            [0.0, 0.0, 0.0], [10.0, 10.0, 10.0], periodic=(True, True, True)
        )
        assert "periodic=(true, true, true)" in repr(c)

    def test_inside_sphere(self):
        c = molpack.InsideSphereRestraint([0.0, 0.0, 0.0], 5.0)
        assert "InsideSphereRestraint" in repr(c)

    def test_outside_sphere(self):
        c = molpack.OutsideSphereRestraint([0.0, 0.0, 0.0], 2.0)
        assert "OutsideSphereRestraint" in repr(c)

    def test_above_plane(self):
        c = molpack.AbovePlaneRestraint([0.0, 0.0, 1.0], 0.0)
        assert "AbovePlaneRestraint" in repr(c)

    def test_below_plane(self):
        c = molpack.BelowPlaneRestraint([0.0, 0.0, 1.0], 10.0)
        assert "BelowPlaneRestraint" in repr(c)


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
            molpack.Target(frame, count=1)
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
            )
            .with_restraint(molpack.InsideSphereRestraint([5.0, 5.0, 5.0], 5.0))
        )
        assert t is not None

    def test_three_restraints(self):
        frame = self._make_frame()
        t = (
            molpack.Target(frame, count=1)
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
            )
            .with_restraint(molpack.InsideSphereRestraint([5.0, 5.0, 5.0], 5.0))
            .with_restraint(molpack.AbovePlaneRestraint([0.0, 0.0, 1.0], 0.0))
        )
        assert t is not None

    def test_all_restraint_types(self):
        frame = self._make_frame()
        t = molpack.Target(frame, count=1)
        for r in [
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
            molpack.InsideSphereRestraint([0.0, 0.0, 0.0], 5.0),
            molpack.OutsideSphereRestraint([0.0, 0.0, 0.0], 2.0),
            molpack.AbovePlaneRestraint([0.0, 0.0, 1.0], 0.0),
            molpack.BelowPlaneRestraint([0.0, 0.0, 1.0], 10.0),
        ]:
            t = t.with_restraint(r)
        assert t is not None
