import numpy as np
import pytest

import molpack


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
        t = molpack.Target("mol", _make_frame(2), 10)
        assert t.natoms == 2
        assert t.count == 10

    def test_elements_from_frame(self):
        t = molpack.Target("water", _make_frame(3, ["O", "H", "H"]), 5)
        assert t.elements == ["O", "H", "H"]

    def test_is_fixed_default_false(self):
        t = molpack.Target("mol", _make_frame(), 1)
        assert t.is_fixed is False

    def test_missing_atoms_block_raises(self):
        with pytest.raises((ValueError, KeyError)):
            molpack.Target("mol", {}, 1)

    def test_missing_x_column_raises(self):
        frame = {
            "atoms": {
                "y": np.array([0.0]),
                "z": np.array([0.0]),
                "element": ["X"],
            }
        }
        with pytest.raises(ValueError, match='"x"'):
            molpack.Target("mol", frame, 1)

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
            molpack.Target("mol", frame, 1)

    def test_accepts_plain_dict_frame(self):
        frame = {
            "atoms": {
                "x": [0.0, 1.0],
                "y": [0.0, 0.0],
                "z": [0.0, 0.0],
                "element": ["C", "H"],
            }
        }
        t = molpack.Target("mol", frame, 3)
        assert t.natoms == 2
        assert t.elements == ["C", "H"]


class TestTargetBuilder:
    def _make_target(self) -> molpack.Target:
        return molpack.Target("mol", _make_frame(), 5)

    def _make_two_atom_target(self) -> molpack.Target:
        return molpack.Target("mol", _make_two_atom_frame(), 5)

    def test_with_name_immutable(self):
        t = self._make_target()
        t2 = t.with_name("water")
        assert "water" in repr(t2)
        assert "water" not in repr(t)

    def test_with_restraint(self):
        t = self._make_target()
        c = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        t2 = t.with_restraint(c)
        assert t2 is not t

    def test_with_multiple_restraints(self):
        t = self._make_target()
        t2 = t.with_restraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        ).with_restraint(molpack.OutsideSphere(2.0, [10.0, 10.0, 10.0]))
        assert t2 is not t

    def test_with_restraint_type_error(self):
        t = self._make_target()
        with pytest.raises(TypeError):
            t.with_restraint("not_a_constraint")  # ty: ignore[invalid-argument-type]

    def test_with_restraint_for_atoms(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        t2 = t.with_restraint_for_atoms([1], c)
        assert t2 is not t

    def test_with_restraint_for_atoms_validates_packmol_indices(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        with pytest.raises(ValueError, match="1-based indexing"):
            t.with_restraint_for_atoms([0], c)
        with pytest.raises(ValueError, match="1-based indexing"):
            t.with_restraint_for_atoms([3], c)

    def test_with_maxmove(self):
        t = self._make_target()
        t2 = t.with_maxmove(100)
        assert t2 is not t

    def test_with_center(self):
        t = self._make_target()
        t2 = t.with_center()
        assert t2 is not t

    def test_without_centering(self):
        t = self._make_target()
        t2 = t.without_centering()
        assert t2 is not t

    def test_constrain_rotation(self):
        t = self._make_target()
        tx = t.constrain_rotation_x(0.0, 10.0)
        ty = t.constrain_rotation_y(30.0, 5.0)
        tz = t.constrain_rotation_z(90.0, 15.0)
        assert tx is not t
        assert ty is not t
        assert tz is not t

    def test_fixed_at(self):
        t = self._make_target()
        t2 = t.fixed_at([5.0, 5.0, 5.0])
        assert t2.count == 1
        assert t2.is_fixed is True
        assert t.is_fixed is False

    def test_fixed_at_with_euler(self):
        t = self._make_target()
        t2 = t.fixed_at_with_euler([5.0, 5.0, 5.0], [45.0, 0.0, 90.0])
        assert t2.count == 1
        assert t2.is_fixed is True
        assert t.is_fixed is False

    def test_repr(self):
        t = molpack.Target("test_mol", _make_frame(), 5)
        r = repr(t)
        assert "Target" in r
        assert "test_mol" in r
        assert "natoms=1" in r


class TestMolpack:
    def test_creation_defaults(self):
        p = molpack.Molpack()
        r = repr(p)
        assert "2.00" in r
        assert "0.0100" in r

    def test_creation_custom(self):
        p = molpack.Molpack(tolerance=3.0, precision=0.1, nloop0=50, sidemax=500.0)
        r = repr(p)
        assert "3.00" in r
        assert "0.1000" in r
        assert "50" in r
        assert "500.0" in r

    def test_builder_immutability(self):
        p1 = molpack.Molpack(tolerance=2.0)
        p2 = p1.with_tolerance(3.0)
        p3 = p1.with_precision(0.5)
        p4 = p1.with_maxit(50)
        p5 = p1.with_nloop0(40)
        p6 = p1.with_sidemax(200.0)
        p7 = p1.with_movefrac(0.1)
        p8 = p1.with_movebadrandom(True)
        p9 = p1.with_disable_movebad(True)
        p10 = p1.with_pbc([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        p11 = p1.with_pbc_box([20.0, 20.0, 20.0])
        p12 = p1.with_progress(False)
        assert "2.00" in repr(p1)
        assert "3.00" in repr(p2)
        assert "0.5000" in repr(p3)
        assert p4 is not p1
        assert p5 is not p1
        assert p6 is not p1
        assert p7 is not p1
        assert p8 is not p1
        assert p9 is not p1
        assert p10 is not p1
        assert p11 is not p1
        assert p12 is not p1


class TestMolpackPack:
    def _make_target(self, count: int = 3) -> molpack.Target:
        return molpack.Target("mol", _make_frame(), count).with_restraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

    def test_minimal_packing(self):
        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        result = packer.pack([self._make_target()], max_loops=50, seed=42)

        assert result.positions.shape == (3, 3)
        assert result.fdist >= 0.0
        assert result.frest >= 0.0
        assert isinstance(result.converged, bool)

    def test_result_elements_match_positions(self):
        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        result = packer.pack([self._make_target()], max_loops=50, seed=42)

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
        target = molpack.Target("water", frame, count=2).with_restraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=50, seed=42)

        assert len(result.elements) == 6
        assert result.positions.shape[0] == 6
        assert "X" not in result.elements
        assert set(result.elements) == {"O", "H"}

    def test_result_elements_order_multiple_targets(self):
        box_c = molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
        t1 = molpack.Target("a", _make_frame(1), 2).with_restraint(box_c)
        t2 = molpack.Target("b", _make_frame(2, ["C", "H"]), 3).with_restraint(box_c)

        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        result = packer.pack([t1, t2], max_loops=50, seed=42)

        assert len(result.elements) == 8
        assert result.positions.shape[0] == 8

    def test_with_seed_reproducible(self):
        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        r1 = packer.pack([self._make_target()], max_loops=30, seed=123)
        r2 = packer.pack([self._make_target()], max_loops=30, seed=123)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_no_targets_raises_runtime_error(self):
        packer = molpack.Molpack().with_progress(False)
        with pytest.raises(RuntimeError, match="No targets"):
            packer.pack([], max_loops=10)

    def test_multiple_targets(self):
        box_constraint = molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        t1 = molpack.Target("a", _make_frame(), 2).with_restraint(box_constraint)
        t2 = molpack.Target("b", _make_frame(), 3).with_restraint(box_constraint)

        packer = molpack.Molpack(tolerance=2.0).with_progress(False)
        result = packer.pack([t1, t2], max_loops=50, seed=42)
        assert result.positions.shape[0] == 5

    def test_pack_result_repr(self):
        packer = molpack.Molpack().with_progress(False)
        result = packer.pack([self._make_target(2)], max_loops=30, seed=1)
        r = repr(result)
        assert "PackResult" in r
        assert "converged" in r

    def test_extended_target_and_packer_options(self):
        frame = {
            "atoms": {
                "x": np.array([0.0, 1.5], dtype=np.float64),
                "y": np.zeros(2, dtype=np.float64),
                "z": np.zeros(2, dtype=np.float64),
                "element": ["C", "H"],
            }
        }
        target = (
            molpack.Target("mol", frame, 2)
            .with_restraint(molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]))
            .with_restraint_for_atoms([1], molpack.AbovePlane([0.0, 0.0, 1.0], 0.0))
            .constrain_rotation_x(0.0, 30.0)
            .constrain_rotation_y(0.0, 30.0)
            .constrain_rotation_z(0.0, 30.0)
        )

        packer = molpack.Molpack(
            tolerance=2.0,
            precision=0.01,
            nloop0=20,
            sidemax=100.0,
            movefrac=0.05,
            movebadrandom=True,
            disable_movebad=False,
            progress=False,
        ).with_pbc_box([20.0, 20.0, 20.0])
        result = packer.pack([target], max_loops=30, seed=42)

        assert result.positions.shape == (4, 3)
        assert result.natoms == 4
