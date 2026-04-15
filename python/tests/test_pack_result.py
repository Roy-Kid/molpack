"""Unit tests for ``PackResult`` getters and invariants."""

from __future__ import annotations

import numpy as np
import pytest

import molpack


def _make_tiny_pack():
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    radii = np.array([1.0, 1.0], dtype=np.float64)
    target = molpack.Target.from_coords(
        positions, radii, 3, elements=["O", "H"]
    ).with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0]))
    packer = molpack.Packer(tolerance=2.0).with_progress(False)
    return packer.pack([target], max_loops=50, seed=42)


class TestPackResultProperties:
    def test_positions_dtype_is_float64(self):
        result = _make_tiny_pack()
        assert result.positions.dtype == np.float64

    def test_positions_shape(self):
        result = _make_tiny_pack()
        # 3 copies × 2 atoms = 6
        assert result.positions.shape == (6, 3)

    def test_elements_is_list_of_str(self):
        result = _make_tiny_pack()
        assert isinstance(result.elements, list)
        assert all(isinstance(e, str) for e in result.elements)

    def test_natoms_matches_shape(self):
        result = _make_tiny_pack()
        assert result.natoms == result.positions.shape[0]
        assert result.natoms == len(result.elements)

    def test_converged_is_bool(self):
        result = _make_tiny_pack()
        assert isinstance(result.converged, bool)

    def test_fdist_frest_nonneg(self):
        result = _make_tiny_pack()
        assert result.fdist >= 0.0
        assert result.frest >= 0.0

    def test_element_pattern_repeats_per_copy(self):
        result = _make_tiny_pack()
        # Template is ["O", "H"] × 3 copies → "OHOHOH".
        assert result.elements == ["O", "H"] * 3


class TestPackerErrorPaths:
    def test_empty_targets_list_raises(self):
        packer = molpack.Packer().with_progress(False)
        with pytest.raises(RuntimeError):
            packer.pack([], max_loops=10, seed=1)

    def test_invalid_pbc_raises_runtime(self):
        # Zero-length axis is rejected inside pack().
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        target = molpack.Target.from_coords(positions, radii, 1).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        )
        packer = molpack.Packer().with_progress(False).with_pbc_box([0.0, 10.0, 10.0])
        with pytest.raises(RuntimeError):
            packer.pack([target], max_loops=10, seed=1)


class TestPackerAcceptsCompositeConstraint:
    def test_molecule_constraint_bundle(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([1.0], dtype=np.float64)
        bundle = molpack.InsideBox([0.0, 0.0, 0.0], [20.0, 20.0, 20.0]).and_(
            molpack.OutsideSphere(2.0, [10.0, 10.0, 10.0])
        )

        target = molpack.Target.from_coords(positions, radii, 3).with_constraint(bundle)
        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=100, seed=42)

        assert result.positions.shape == (3, 3)
