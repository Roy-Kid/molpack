"""Unit tests for ``State`` getters and invariants."""

from __future__ import annotations

import math

import molrs
import numpy as np

import molpack


def _col(frame, block: str, name: str) -> np.ndarray:
    """Read a column from a ``molrs.Frame`` block as a numpy array."""
    return np.asarray(frame[block][name])


def _make_frame(
    positions: np.ndarray,
    elements: list[str],
) -> molrs.Frame:
    return molrs.Frame(
        {
            "atoms": {
                "x": positions[:, 0].copy(),
                "y": positions[:, 1].copy(),
                "z": positions[:, 2].copy(),
                "element": elements,
            }
        }
    )


def _make_tiny_pack() -> molpack.State:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    frame = _make_frame(positions, ["O", "H"])
    target = molpack.Target(frame, 3).with_restraint(
        molrs.Cuboid([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
    )
    packer = molpack.GenCanPack().with_tolerance(2.0).with_progress(False).with_seed(42)
    return packer.run([target], max_loops=50)


class TestState:
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
        assert isinstance(frame, molrs.Frame)
        for col in ("x", "y", "z", "element", "id", "mol_id"):
            assert len(_col(frame, "atoms", col)) == result.natoms

    def test_repr_starts_with_state(self):
        result = _make_tiny_pack()
        assert repr(result).startswith("State(")


class TestIntraResidual:
    """``State.intra`` forwards nested scored/exempted (Å); never aliases."""

    def test_intra_is_intra_residual_with_scored_exempted_floats(self):
        result = _make_tiny_pack()
        intra = result.intra
        assert type(intra).__name__ == "IntraResidual"
        assert isinstance(intra, molpack.IntraResidual)
        assert isinstance(intra.scored, float)
        assert isinstance(intra.exempted, float)

    def test_no_min_intra_aliases(self):
        result = _make_tiny_pack()
        assert not hasattr(result, "min_intra_scored")
        assert not hasattr(result, "min_intra_exempt")

    def test_bonded_diatomic_scored_is_infinite(self):
        result = TestFrameTopology()._pack(1)
        assert result.intra.scored == math.inf


class TestFrameTopology:
    """End-to-end: a template's topology is replayed onto packed coordinates."""

    @staticmethod
    def _diatomic_with_bond() -> molrs.Frame:
        return molrs.Frame(
            {
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
        )

    def _pack(self, copies: int, box: bool = False) -> molpack.State:
        target = molpack.Target(self._diatomic_with_bond(), copies).with_restraint(
            molrs.Cuboid([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        )
        packer = (
            molpack.GenCanPack().with_tolerance(2.0).with_progress(False).with_seed(7)
        )
        if box:
            packer = packer.with_periodic_box([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        return packer.run([target], max_loops=50)

    def test_frame_carries_replicated_topology(self):
        result = self._pack(3)
        frame = result.frame

        assert isinstance(frame, molrs.Frame)
        assert np.array_equal(_col(frame, "atoms", "id"), np.arange(1, 7))
        assert np.array_equal(
            _col(frame, "atoms", "mol_id"), np.array([1, 1, 2, 2, 3, 3])
        )
        assert np.array_equal(_col(frame, "atoms", "type"), np.array(["A", "B"] * 3))

    def test_bond_indices_offset_per_copy(self):
        frame = self._pack(3).frame
        assert np.array_equal(_col(frame, "bonds", "atomi"), np.array([0, 2, 4]))
        assert np.array_equal(_col(frame, "bonds", "atomj"), np.array([1, 3, 5]))
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

    def test_frame_getter_returns_same_object(self):
        result = self._pack(3)
        assert result.frame is result.frame

    def test_assigned_box_persists_on_frame(self):
        result = self._pack(3)
        result.frame.box = molrs.Box.cube(20.0)
        assert result.frame.box is not None
        assert np.allclose(np.asarray(result.frame.box.lengths), [20.0, 20.0, 20.0])
