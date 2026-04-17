"""Unit tests for ``PackResult`` getters and invariants."""

from __future__ import annotations

import numpy as np
import pytest

import molpack


def _make_frame(
    positions: np.ndarray,
    elements: list[str],
) -> dict:
    return {
        "atoms": {
            "x": positions[:, 0].copy(),
            "y": positions[:, 1].copy(),
            "z": positions[:, 2].copy(),
            "element": elements,
        }
    }


def _make_tiny_pack() -> molpack.PackResult:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    frame = _make_frame(positions, ["O", "H"])
    target = molpack.Target(frame, 3).with_restraint(
        molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
    )
    packer = molpack.Molpack().with_tolerance(2.0).with_progress(False).with_seed(42)
    return packer.pack([target], max_loops=50)


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

    def test_frame_exposes_atoms_block(self):
        result = _make_tiny_pack()
        frame = result.frame
        assert "atoms" in frame
        atoms = frame["atoms"]
        assert "x" in atoms and "y" in atoms and "z" in atoms
        assert "element" in atoms
        assert len(atoms["element"]) == result.natoms
        assert atoms["x"].shape[0] == result.natoms


class TestMolpackErrorPaths:
    def test_empty_targets_list_raises(self):
        packer = molpack.Molpack().with_progress(False).with_seed(1)
        with pytest.raises(molpack.NoTargetsError):
            packer.pack([], max_loops=10)

    def test_invalid_pbc_raises_typed_error(self):
        # Zero-length axis on a periodic box is rejected at pack().
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        frame = {
            "atoms": {
                "x": positions[:, 0],
                "y": positions[:, 1],
                "z": positions[:, 2],
                "element": ["X"],
            }
        }
        target = molpack.Target(frame, 1).with_restraint(
            molpack.InsideBoxRestraint(
                [0.0, 0.0, 0.0], [0.0, 10.0, 10.0], periodic=(True, True, True)
            )
        )
        packer = molpack.Molpack().with_progress(False).with_seed(1)
        with pytest.raises(molpack.InvalidPBCBoxError):
            packer.pack([target], max_loops=10)

    def test_pack_error_is_runtime_error_subclass(self):
        assert issubclass(molpack.NoTargetsError, molpack.PackError)
        assert issubclass(molpack.PackError, RuntimeError)


class TestMultipleRestraints:
    def test_stacked_restraints(self):
        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        frame = {
            "atoms": {
                "x": positions[:, 0],
                "y": positions[:, 1],
                "z": positions[:, 2],
                "element": ["X"],
            }
        }
        target = (
            molpack.Target(frame, 3)
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
            )
            .with_restraint(molpack.OutsideSphereRestraint([10.0, 10.0, 10.0], 2.0))
        )
        packer = (
            molpack.Molpack().with_tolerance(2.0).with_progress(False).with_seed(42)
        )
        result = packer.pack([target], max_loops=100)

        assert result.positions.shape == (3, 3)
