import numpy as np
import pytest

import molpack
from molpack import Angle, Axis, CenteringMode


def _make_frame(
    n: int = 1,
    elements: list[str] | None = None,
) -> dict:
    """Minimal frame dict with n atoms at the origin."""
    if elements is None:
        elements = ["X"] * n
    positions = np.zeros((n, 3), dtype=np.float64)
    return {
        "atoms": {
            "x": positions[:, 0].copy(),
            "y": positions[:, 1].copy(),
            "z": positions[:, 2].copy(),
            "element": elements,
        }
    }


def _make_two_atom_frame() -> dict:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    return {
        "atoms": {
            "x": positions[:, 0].copy(),
            "y": positions[:, 1].copy(),
            "z": positions[:, 2].copy(),
            "element": ["O", "H"],
        }
    }


class TestTargetConstructor:
    def test_basic(self):
        t = molpack.Target(_make_frame(2), 10)
        assert t.natoms == 2
        assert t.count == 10
        assert t.name is None

    def test_with_name(self):
        t = molpack.Target(_make_frame(2), 10).with_name("mol")
        assert t.name == "mol"

    def test_elements_from_frame(self):
        t = molpack.Target(_make_frame(3, ["O", "H", "H"]), 5)
        assert t.elements == ["O", "H", "H"]

    def test_is_fixed_default_false(self):
        t = molpack.Target(_make_frame(), 1)
        assert t.is_fixed is False

    def test_missing_atoms_block_raises(self):
        with pytest.raises((ValueError, KeyError)):
            molpack.Target({}, 1)

    def test_missing_x_column_raises(self):
        frame = {
            "atoms": {
                "y": np.array([0.0]),
                "z": np.array([0.0]),
                "element": ["X"],
            }
        }
        with pytest.raises(ValueError, match='"x"'):
            molpack.Target(frame, 1)

    def test_mismatched_lengths_raises(self):
        frame = {
            "atoms": {
                "x": np.array([0.0, 1.0], dtype=np.float64),
                "y": np.array([0.0], dtype=np.float64),
                "z": np.array([0.0, 0.0], dtype=np.float64),
                "element": ["X", "X"],
            }
        }
        with pytest.raises(ValueError):
            molpack.Target(frame, 1)

    def test_accepts_plain_dict_frame(self):
        frame = {
            "atoms": {
                "x": [0.0, 1.0],
                "y": [0.0, 0.0],
                "z": [0.0, 0.0],
                "element": ["C", "H"],
            }
        }
        t = molpack.Target(frame, 3)
        assert t.natoms == 2
        assert t.elements == ["C", "H"]


class TestTargetBuilder:
    def _make_target(self) -> molpack.Target:
        return molpack.Target(_make_frame(), 5)

    def _make_two_atom_target(self) -> molpack.Target:
        return molpack.Target(_make_two_atom_frame(), 5)

    def test_with_name_immutable(self):
        t = self._make_target()
        t2 = t.with_name("water")
        assert t2.name == "water"
        assert t.name is None

    def test_with_restraint(self):
        t = self._make_target()
        c = molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        t2 = t.with_restraint(c)
        assert t2 is not t

    def test_with_multiple_restraints(self):
        t = self._make_target()
        t2 = t.with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        ).with_restraint(molpack.OutsideSphereRestraint([10.0, 10.0, 10.0], 2.0))
        assert t2 is not t

    def test_with_restraint_type_error(self):
        t = self._make_target()
        with pytest.raises(TypeError):
            t.with_restraint("not_a_constraint")  # ty: ignore[invalid-argument-type]

    def test_with_atom_restraint(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphereRestraint([0.0, 0.0, 0.0], 5.0)
        # 0-based: atom index 0 is the first atom.
        t2 = t.with_atom_restraint([0], c)
        assert t2 is not t

    def test_with_atom_restraint_rejects_out_of_range(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphereRestraint([0.0, 0.0, 0.0], 5.0)
        with pytest.raises(ValueError, match="0-based"):
            t.with_atom_restraint([2], c)

    def test_with_perturb_budget(self):
        t = self._make_target()
        t2 = t.with_perturb_budget(100)
        assert t2 is not t

    def test_with_centering_center(self):
        t = self._make_target()
        t2 = t.with_centering(CenteringMode.CENTER)
        assert t2 is not t

    def test_with_centering_off(self):
        t = self._make_target()
        t2 = t.with_centering(CenteringMode.OFF)
        assert t2 is not t

    def test_with_rotation_bound(self):
        t = self._make_target()
        tx = t.with_rotation_bound(
            Axis.X, Angle.from_degrees(0.0), Angle.from_degrees(10.0)
        )
        ty = t.with_rotation_bound(
            Axis.Y, Angle.from_degrees(30.0), Angle.from_degrees(5.0)
        )
        tz = t.with_rotation_bound(
            Axis.Z, Angle.from_degrees(90.0), Angle.from_degrees(15.0)
        )
        assert tx is not t
        assert ty is not t
        assert tz is not t

    def test_fixed_at(self):
        t = molpack.Target(_make_frame(), 1)
        t2 = t.fixed_at([5.0, 5.0, 5.0])
        assert t2.count == 1
        assert t2.is_fixed is True
        assert t.is_fixed is False

    def test_fixed_at_with_orientation(self):
        t = molpack.Target(_make_frame(), 1)
        t2 = t.fixed_at([5.0, 5.0, 5.0]).with_orientation(
            (
                Angle.from_degrees(45.0),
                Angle.ZERO,
                Angle.from_degrees(90.0),
            )
        )
        assert t2.count == 1
        assert t2.is_fixed is True

    def test_repr(self):
        t = molpack.Target(_make_frame(), 5).with_name("test_mol")
        r = repr(t)
        assert "Target" in r
        assert "test_mol" in r
        assert "natoms=1" in r


class TestMolpack:
    def test_creation_default(self):
        p = molpack.Molpack()
        r = repr(p)
        assert "Molpack" in r

    def test_builder_immutability(self):
        p1 = molpack.Molpack()
        p2 = p1.with_tolerance(3.0)
        p3 = p1.with_precision(0.5)
        p4 = p1.with_inner_iterations(50)
        p5 = p1.with_init_passes(40)
        p6 = p1.with_init_box_half_size(200.0)
        p7 = p1.with_perturb_fraction(0.1)
        p8 = p1.with_random_perturb(True)
        p9 = p1.with_perturb(False)
        p10 = p1.with_progress(False)
        p11 = p1.with_seed(42)
        p12 = p1.with_parallel_eval(True)
        for later in (p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12):
            assert later is not p1


class TestMolpackPack:
    def _make_target(self, count: int = 3) -> molpack.Target:
        return molpack.Target(_make_frame(), count).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

    def _packer(self) -> molpack.Molpack:
        return molpack.Molpack().with_tolerance(2.0).with_progress(False)

    def test_minimal_packing(self):
        result = self._packer().with_seed(42).pack([self._make_target()], max_loops=50)

        assert result.positions.shape == (3, 3)
        assert result.fdist >= 0.0
        assert result.frest >= 0.0
        assert isinstance(result.converged, bool)

    def test_result_elements_match_positions(self):
        result = self._packer().with_seed(42).pack([self._make_target()], max_loops=50)

        assert len(result.elements) == result.positions.shape[0]
        assert result.natoms == result.positions.shape[0]
        assert all(e == "X" for e in result.elements)

    def test_result_elements_from_frame(self):
        frame = {
            "atoms": {
                "x": np.array([0.0, 0.96, -0.24], dtype=np.float64),
                "y": np.array([0.0, 0.0, 0.93], dtype=np.float64),
                "z": np.zeros(3, dtype=np.float64),
                "element": ["O", "H", "H"],
            }
        }
        target = (
            molpack.Target(frame, count=2)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
            )
        )

        result = self._packer().with_seed(42).pack([target], max_loops=50)

        assert len(result.elements) == 6
        assert result.positions.shape[0] == 6
        assert "X" not in result.elements
        assert set(result.elements) == {"O", "H"}

    def test_result_elements_order_multiple_targets(self):
        box_c = molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
        t1 = molpack.Target(_make_frame(1), 2).with_restraint(box_c)
        t2 = molpack.Target(_make_frame(2, ["C", "H"]), 3).with_restraint(box_c)

        result = self._packer().with_seed(42).pack([t1, t2], max_loops=50)

        assert len(result.elements) == 8
        assert result.positions.shape[0] == 8

    def test_with_seed_reproducible(self):
        packer = self._packer().with_seed(123)
        r1 = packer.pack([self._make_target()], max_loops=30)
        r2 = packer.pack([self._make_target()], max_loops=30)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_no_targets_raises_typed_error(self):
        packer = molpack.Molpack().with_progress(False)
        with pytest.raises(molpack.NoTargetsError):
            packer.pack([], max_loops=10)

    def test_multiple_targets(self):
        box = molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        t1 = molpack.Target(_make_frame(), 2).with_restraint(box)
        t2 = molpack.Target(_make_frame(), 3).with_restraint(box)

        result = self._packer().with_seed(42).pack([t1, t2], max_loops=50)
        assert result.positions.shape[0] == 5

    def test_pack_result_repr(self):
        result = self._packer().with_seed(1).pack([self._make_target(2)], max_loops=30)
        r = repr(result)
        assert "PackResult" in r
        assert "converged" in r


class TestMolpackGlobalRestraint:
    def test_global_restraint_broadcasts(self):
        frame = _make_frame()
        t1 = molpack.Target(frame, 2).with_name("a")
        t2 = molpack.Target(frame, 2).with_name("b")

        packer = (
            molpack.Molpack()
            .with_progress(False)
            .with_global_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
            )
            .with_seed(42)
        )
        result = packer.pack([t1, t2], max_loops=50)
        assert result.natoms == 4
        # All atoms must land inside the global box.
        for pos in result.positions:
            for k in range(3):
                assert -0.1 <= pos[k] <= 30.1
