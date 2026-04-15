import numpy as np
import pytest

import molpack


class TestTargetFromCoords:
    def test_basic(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0, 1.0], dtype=np.float64)
        t = molpack.Target.from_coords(positions, radii, 10)
        assert t.natoms == 2
        assert t.count == 10

    def test_elements_default_to_x(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        t = molpack.Target.from_coords(positions, radii, 1)
        assert t.elements == ["X"]

    def test_elements_override(self):
        positions = np.array(
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
            dtype=np.float64,
        )
        radii = np.array([1.52, 1.20, 1.20], dtype=np.float64)
        t = molpack.Target.from_coords(positions, radii, 5, elements=["O", "H", "H"])
        assert t.elements == ["O", "H", "H"]

    def test_is_fixed_default_false(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        t = molpack.Target.from_coords(positions, radii, 1)
        assert t.is_fixed is False

    def test_bad_positions_shape(self):
        with pytest.raises(ValueError, match="N, 3"):
            molpack.Target.from_coords(
                np.ones((3, 2), dtype=np.float64),
                np.ones(3, dtype=np.float64),
                5,
            )

    def test_mismatched_lengths(self):
        with pytest.raises(ValueError, match="same length"):
            molpack.Target.from_coords(
                np.ones((3, 3), dtype=np.float64),
                np.ones(2, dtype=np.float64),
                5,
            )

    def test_elements_length_mismatch(self):
        with pytest.raises(ValueError, match="elements"):
            molpack.Target.from_coords(
                np.ones((2, 3), dtype=np.float64),
                np.ones(2, dtype=np.float64),
                1,
                elements=["O"],
            )


class TestTargetBuilder:
    def _make_target(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        return molpack.Target.from_coords(positions, radii, 5)

    def _make_two_atom_target(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0, 1.0], dtype=np.float64)
        return molpack.Target.from_coords(positions, radii, 5)

    def test_with_name_immutable(self):
        t = self._make_target()
        t2 = t.with_name("water")
        assert "water" in repr(t2)
        assert "water" not in repr(t)

    def test_with_constraint(self):
        t = self._make_target()
        c = molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        t2 = t.with_constraint(c)
        assert t2 is not t

    def test_with_constraint_type_error(self):
        t = self._make_target()
        with pytest.raises(TypeError):
            t.with_constraint("not_a_constraint")

    def test_with_constraint_for_atoms(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        t2 = t.with_constraint_for_atoms([1], c)
        assert t2 is not t

    def test_with_constraint_for_atoms_validates_packmol_indices(self):
        t = self._make_two_atom_target()
        c = molpack.InsideSphere(5.0, [0.0, 0.0, 0.0])
        with pytest.raises(ValueError, match="1-based indexing"):
            t.with_constraint_for_atoms([0], c)
        with pytest.raises(ValueError, match="1-based indexing"):
            t.with_constraint_for_atoms([3], c)

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
        t = self._make_target().with_name("test_mol")
        r = repr(t)
        assert "Target" in r
        assert "test_mol" in r
        assert "natoms=1" in r


class TestPacker:
    def test_creation_defaults(self):
        p = molpack.Packer()
        r = repr(p)
        assert "2.00" in r
        assert "0.0100" in r

    def test_creation_custom(self):
        p = molpack.Packer(tolerance=3.0, precision=0.1)
        r = repr(p)
        assert "3.00" in r
        assert "0.1000" in r

    def test_builder_immutability(self):
        p1 = molpack.Packer(tolerance=2.0)
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


class TestPackerPack:
    def test_minimal_packing(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        target = molpack.Target.from_coords(positions, radii, 3).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=50, seed=42)

        assert result.positions.shape == (3, 3)
        assert result.fdist >= 0.0
        assert result.frest >= 0.0
        assert isinstance(result.converged, bool)

    def test_result_elements_match_positions(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        target = molpack.Target.from_coords(positions, radii, 3).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=50, seed=42)

        assert len(result.elements) == result.positions.shape[0]
        assert result.natoms == result.positions.shape[0]
        assert all(e == "X" for e in result.elements)

    def test_result_elements_from_override(self):
        positions = np.array(
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
            dtype=np.float64,
        )
        radii = np.array([1.52, 1.20, 1.20], dtype=np.float64)
        target = molpack.Target.from_coords(
            positions, radii, 2, elements=["O", "H", "H"]
        ).with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]))

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=50, seed=42)

        # 2 water molecules * 3 atoms = 6 atoms
        assert len(result.elements) == 6
        assert result.positions.shape[0] == 6
        assert "X" not in result.elements
        assert set(result.elements) == {"O", "H"}

    def test_result_elements_order_multiple_targets(self):
        pos1 = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        rad1 = np.array([1.0], dtype=np.float64)
        pos2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
        rad2 = np.array([1.0, 1.0], dtype=np.float64)

        box_c = molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
        t1 = molpack.Target.from_coords(pos1, rad1, 2).with_constraint(box_c)
        t2 = molpack.Target.from_coords(pos2, rad2, 3).with_constraint(box_c)

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([t1, t2], max_loops=50, seed=42)

        assert len(result.elements) == 8
        assert result.positions.shape[0] == 8

    def test_with_seed_reproducible(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        target = molpack.Target.from_coords(positions, radii, 3).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        r1 = packer.pack([target], max_loops=30, seed=123)
        r2 = packer.pack([target], max_loops=30, seed=123)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_no_targets_raises_runtime_error(self):
        packer = molpack.Packer().with_progress(False)
        with pytest.raises(RuntimeError, match="No targets"):
            packer.pack([], max_loops=10)

    def test_multiple_targets(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        box_constraint = molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])

        t1 = molpack.Target.from_coords(positions, radii, 2).with_constraint(
            box_constraint
        )
        t2 = molpack.Target.from_coords(positions, radii, 3).with_constraint(
            box_constraint
        )

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([t1, t2], max_loops=50, seed=42)
        assert result.positions.shape[0] == 5

    def test_pack_result_repr(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        target = molpack.Target.from_coords(positions, radii, 2).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )
        packer = molpack.Packer().with_progress(False)
        result = packer.pack([target], max_loops=30, seed=1)
        r = repr(result)
        assert "PackResult" in r
        assert "converged" in r

    def test_extended_target_and_packer_options(self):
        positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0, 1.0], dtype=np.float64)
        target = (
            molpack.Target.from_coords(positions, radii, 2)
            .with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]))
            .with_constraint_for_atoms([1], molpack.AbovePlane([0.0, 0.0, 1.0], 0.0))
            .constrain_rotation_x(0.0, 30.0)
            .constrain_rotation_y(0.0, 30.0)
            .constrain_rotation_z(0.0, 30.0)
        )

        packer = (
            molpack.Packer(tolerance=2.0, precision=0.01)
            .with_nloop0(20)
            .with_sidemax(100.0)
            .with_movefrac(0.05)
            .with_movebadrandom(True)
            .with_disable_movebad(False)
            .with_pbc_box([20.0, 20.0, 20.0])
            .with_progress(False)
        )
        result = packer.pack([target], max_loops=30, seed=42)

        assert result.positions.shape == (4, 3)
        assert result.natoms == 4
