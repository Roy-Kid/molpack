import pytest

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


class TestConstraintComposition:
    def test_two_constraints(self):
        c1 = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        c2 = molpack.InsideSphere(5.0, [5.0, 5.0, 5.0])
        combined = c1.and_(c2)
        assert "restraints=2" in repr(combined)

    def test_chain_three(self):
        c1 = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        c2 = molpack.InsideSphere(5.0, [5.0, 5.0, 5.0])
        c3 = molpack.AbovePlane([0.0, 0.0, 1.0], 0.0)
        combined = c1.and_(c2).and_(c3)
        assert "restraints=3" in repr(combined)

    def test_and_each_type(self):
        box_ = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        sphere = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        outside = molpack.OutsideSphere(2.0, [0.0, 0.0, 0.0])
        above = molpack.AbovePlane([0.0, 0.0, 1.0], 0.0)
        below = molpack.BelowPlane([0.0, 0.0, 1.0], 10.0)

        assert box_.and_(sphere) is not None
        assert sphere.and_(outside) is not None
        assert outside.and_(above) is not None
        assert above.and_(below) is not None
        assert below.and_(box_) is not None

    def test_type_error_on_bad_argument(self):
        c = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        with pytest.raises(TypeError):
            c.and_("not_a_constraint")  # ty: ignore[invalid-argument-type]

    def test_molecule_constraint_and(self):
        c1 = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        c2 = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        mc = c1.and_(c2)
        c3 = molpack.AbovePlane([0.0, 0.0, 1.0], 0.0)
        mc2 = mc.and_(c3)
        assert "restraints=3" in repr(mc2)
