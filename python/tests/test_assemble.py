"""Unit tests for ``molpack._assemble`` — the frame assembly module.

These exercise the topology-replay logic directly (independent of a real pack
run) so the per-copy index offsetting and the size-mismatch guard are pinned.
"""

from __future__ import annotations

import molrs
import numpy as np
import pytest

from molpack import _assemble


def _col(frame, block: str, name: str) -> np.ndarray:
    return np.asarray(frame[block].view(name))


def _diatomic() -> dict:
    """A 2-atom template with one 0-based bond and a carried ``type`` column."""
    return {
        "atoms": {
            "type": np.array(["A", "B"]),
            "element": np.array(["C", "H"]),
            "x": np.array([0.0, 1.0]),
            "y": np.zeros(2),
            "z": np.zeros(2),
        },
        "bonds": {"atomi": np.array([0]), "atomj": np.array([1])},
    }


def _monatomic() -> dict:
    return {
        "atoms": {
            "element": np.array(["Ar"]),
            "x": np.zeros(1),
            "y": np.zeros(1),
            "z": np.zeros(1),
        }
    }


class TestTopologyFrame:
    def test_index_columns_offset_per_copy(self):
        positions = np.arange(3 * 2 * 3, dtype=float).reshape(
            -1, 3
        )  # 3 copies × 2 atoms
        frame = _assemble.topology_frame(positions, [(_diatomic(), 3)], None)

        assert isinstance(frame, molrs.Frame)
        assert np.array_equal(_col(frame, "bonds", "atomi"), [0, 2, 4])
        assert np.array_equal(_col(frame, "bonds", "atomj"), [1, 3, 5])
        assert np.array_equal(_col(frame, "bonds", "id"), [1, 2, 3])
        assert np.array_equal(_col(frame, "atoms", "id"), np.arange(1, 7))
        assert np.array_equal(_col(frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 3])
        assert np.array_equal(_col(frame, "atoms", "type"), ["A", "B"] * 3)

    def test_multiple_templates_chain_their_offsets(self):
        # 2 diatomics (4 atoms) then 2 argons (2 atoms) → 6 atoms total.
        positions = np.zeros((6, 3))
        frame = _assemble.topology_frame(
            positions, [(_diatomic(), 2), (_monatomic(), 2)], None
        )
        # Bonds belong only to the diatomics: (0,1) and (2,3).
        assert np.array_equal(_col(frame, "bonds", "atomi"), [0, 2])
        assert np.array_equal(_col(frame, "bonds", "atomj"), [1, 3])
        # mol_id spans all four molecules.
        assert np.array_equal(_col(frame, "atoms", "mol_id"), [1, 1, 2, 2, 3, 4])

    def test_coords_are_taken_from_positions(self):
        positions = np.arange(2 * 2 * 3, dtype=float).reshape(-1, 3)
        frame = _assemble.topology_frame(positions, [(_diatomic(), 2)], None)
        assert np.allclose(_col(frame, "atoms", "x"), positions[:, 0])
        assert np.allclose(_col(frame, "atoms", "z"), positions[:, 2])

    def test_size_mismatch_raises(self):
        positions = np.zeros((5, 3))  # templates need 3×2 = 6
        with pytest.raises(ValueError, match="account for 6"):
            _assemble.topology_frame(positions, [(_diatomic(), 3)], None)

    def test_box_is_stamped(self):
        positions = np.zeros((2, 3))
        box = ([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        frame = _assemble.topology_frame(positions, [(_diatomic(), 1)], box)
        assert frame.box is not None


class TestCoordsOnlyFrame:
    def test_builds_atoms_block(self):
        positions = np.arange(3 * 3, dtype=float).reshape(-1, 3)
        frame = _assemble.coords_only_frame(positions, ["O", "H", "H"], None)
        assert isinstance(frame, molrs.Frame)
        assert np.array_equal(_col(frame, "atoms", "id"), [1, 2, 3])
        assert np.array_equal(_col(frame, "atoms", "element"), ["O", "H", "H"])
        assert np.allclose(_col(frame, "atoms", "y"), positions[:, 1])
