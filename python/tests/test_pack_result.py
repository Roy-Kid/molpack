"""Unit tests for ``PackResult`` getters and invariants."""

from __future__ import annotations

import molrs
import numpy as np
import pytest

import molpack


def _col(frame, block: str, name: str) -> np.ndarray:
    """Read a column from a ``molrs.Frame`` block as a numpy array."""
    return np.asarray(frame[block].view(name))


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

    def test_frame_is_molrs_frame_with_atoms(self):
        result = _make_tiny_pack()
        frame = result.frame
        # `.frame` is a genuine molrs.Frame, not a dict.
        assert isinstance(frame, molrs.Frame)
        for col in ("x", "y", "z", "element", "id", "mol_id"):
            assert len(_col(frame, "atoms", col)) == result.natoms


class TestFrameTopology:
    """End-to-end: a template's topology is replayed onto packed coordinates."""

    @staticmethod
    def _diatomic_with_bond() -> dict:
        return {
            "atoms": {
                "type": np.array(["A", "B"]),
                "charge": np.array([0.1, -0.1]),
                "mass": np.array([12.0, 1.0]),
                "element": np.array(["C", "H"]),
                "x": np.array([0.0, 1.0]),
                "y": np.array([0.0, 0.0]),
                "z": np.array([0.0, 0.0]),
            },
            "bonds": {"atomi": np.array([0]), "atomj": np.array([1])},
        }

    def _pack(self, copies: int, box: bool = False) -> molpack.PackResult:
        target = molpack.Target(self._diatomic_with_bond(), copies).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        )
        packer = molpack.Molpack().with_tolerance(2.0).with_progress(False).with_seed(7)
        if box:
            packer = packer.with_periodic_box([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        return packer.pack([target], max_loops=50)

    def test_frame_carries_replicated_topology(self):
        result = self._pack(3)
        frame = result.frame

        assert isinstance(frame, molrs.Frame)
        # All template columns survive, plus regenerated id / mol_id.
        assert np.array_equal(_col(frame, "atoms", "id"), np.arange(1, 7))
        assert np.array_equal(
            _col(frame, "atoms", "mol_id"), np.array([1, 1, 2, 2, 3, 3])
        )
        assert np.array_equal(_col(frame, "atoms", "type"), np.array(["A", "B"] * 3))

    def test_bond_indices_offset_per_copy(self):
        frame = self._pack(3).frame
        # One bond per copy, second atom of each diatomic: (0,1),(2,3),(4,5).
        assert np.array_equal(_col(frame, "bonds", "atomi"), np.array([0, 2, 4]))
        assert np.array_equal(_col(frame, "bonds", "atomj"), np.array([1, 3, 5]))
        # Regenerated bond ids.
        assert np.array_equal(_col(frame, "bonds", "id"), np.array([1, 2, 3]))

    def test_atom_coords_match_packed_positions(self):
        result = self._pack(3)
        frame = result.frame
        packed = result.positions
        assert np.allclose(_col(frame, "atoms", "x"), packed[:, 0])
        assert np.allclose(_col(frame, "atoms", "y"), packed[:, 1])
        assert np.allclose(_col(frame, "atoms", "z"), packed[:, 2])

    def test_periodic_box_is_stamped_on_frame(self):
        box = self._pack(3, box=True).frame.box
        assert box is not None
        assert np.allclose(np.asarray(box.lengths), [15.0, 15.0, 15.0])

    def test_no_box_when_not_declared(self):
        assert self._pack(3).frame.box is None


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
