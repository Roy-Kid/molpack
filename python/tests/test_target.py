"""``Target.with_special_bonds`` marshalling (special-bonds-06-mirror).

This module never calls ``run()``. Fractional weights store on the Target;
growth refuses them by name in ``test_grow.py``.
"""

from __future__ import annotations

import math

import molrs
import numpy as np
import pytest

from molpack import Target


def _chain_frame(n: int = 5, bond: float = 1.53) -> molrs.store.Frame:
    """Planar zigzag bead chain — local copy of the grow fixture."""
    theta = math.radians(109.5)
    alpha = (math.pi - theta) / 2.0
    dx, dz = bond * math.cos(alpha), bond * math.sin(alpha)
    idx = np.arange(n)
    return molrs.store.Frame(
        {
            "atoms": {
                "x": idx.astype(np.float64) * dx,
                "y": np.zeros(n, dtype=np.float64),
                "z": np.where(idx % 2 == 1, dz, 0.0),
                "element": ["C"] * n,
            },
            "bonds": {
                "atomi": np.arange(0, n - 1, dtype=np.uint64),
                "atomj": np.arange(1, n, dtype=np.uint64),
            },
        }
    )


class TestTarget:
    """Default Cassandra table, builder immutability, and marshalling."""

    def test_default_special_bonds_is_depth_3(self):
        # Cassandra Intra_Scaling depth 3: 1-2/1-3/1-4 exempt, 1-5+ scored.
        t = Target(_chain_frame(), 1)
        assert t.special_bonds == [0.0, 0.0, 0.0, 1.0]

    def test_with_special_bonds_leaves_original_unchanged(self):
        original = Target(_chain_frame(), 1)
        updated = original.with_special_bonds([0.0, 0.0, 1.0])
        assert original.special_bonds == [0.0, 0.0, 0.0, 1.0]
        assert updated is not original
        assert updated.special_bonds == [0.0, 0.0, 1.0]

    def test_cg_table_round_trips(self):
        t = Target(_chain_frame(), 1).with_special_bonds([0.0, 0.0, 1.0])
        assert t.special_bonds == [0.0, 0.0, 1.0]

    def test_empty_list_raises_value_error(self):
        t = Target(_chain_frame(), 1)
        with pytest.raises(ValueError):
            t.with_special_bonds([])

    def test_non_finite_weight_raises_value_error(self):
        t = Target(_chain_frame(), 1)
        with pytest.raises(ValueError):
            t.with_special_bonds([0.0, math.nan, 1.0])
        with pytest.raises(ValueError):
            t.with_special_bonds([0.0, math.inf, 1.0])

    def test_weight_above_one_raises_value_error(self):
        t = Target(_chain_frame(), 1)
        with pytest.raises(ValueError):
            t.with_special_bonds([0.0, 1.5])

    def test_fractional_weight_stores_and_round_trips(self):
        # Amber 1-4 0.5 is a legal table on Target; growth refuses it at run.
        t = Target(_chain_frame(), 1).with_special_bonds([0.0, 0.0, 0.5, 1.0])
        assert t.special_bonds == [0.0, 0.0, 0.5, 1.0]


class TestWithHydrogens:
    """``Target.with_hydrogens`` index marshalling."""

    def test_returns_a_new_target(self):
        original = Target(_chain_frame(), 1)
        assert original.with_hydrogens([0, 4]) is not original

    def test_empty_list_is_accepted(self):
        Target(_chain_frame(), 1).with_hydrogens([])

    def test_out_of_range_index_raises_value_error(self):
        with pytest.raises(ValueError):
            Target(_chain_frame(5), 1).with_hydrogens([5])
