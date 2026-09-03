"""Packs whose only GENCAN-side molecules are fixed (zero movable DOF).

Regression for the `initial.rs` panic (`index out of bounds: the len is 0`)
when the placement vector is empty: a fixed-only pack, and the grow →
fixed-matrix chaining whose rigid stage holds nothing movable.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import molpack

if TYPE_CHECKING:
    import molrs


def _chain3() -> molrs.Frame:
    import molrs

    fr = molrs.Frame()
    fr["atoms"] = {
        "x": np.array([0.0, 1.0, 2.0]),
        "y": np.zeros(3),
        "z": np.zeros(3),
        "element": ["O", "C", "O"],
    }
    fr["bonds"] = {"atomi": np.array([0, 1]), "atomj": np.array([1, 2])}
    return fr


def _box_gencan() -> molpack.GenCanPack:
    return (
        molpack.GenCanPack()
        .with_progress(False)
        .with_seed(7)
        .with_tolerance(2.0)
        .with_periodic_box([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
    )


class TestFixedOnlyPack:
    def test_fixed_only_pack_returns_the_structure(self):
        fixed = (
            molpack.Target(_chain3(), count=1)
            .with_name("matrix")
            .with_centering(molpack.CenteringMode.OFF)
            .fixed_at([0.0, 0.0, 0.0])
        )
        res = _box_gencan().run([fixed], max_loops=5)
        assert res.converged is True
        assert res.natoms == 3

    def test_grown_matrix_fixed_only_rigid_stage(self):
        # The serial-growth composition, explicit form: grow chains, freeze
        # them via ``Target.fixed_from``, and run the rigid stage with the
        # frozen matrix as its ONLY molecule — zero movable DOF again.
        grown = (
            molpack.CbmcGrow(molpack.TorsionPrior.uniform())
            .with_progress(False)
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
            .run([molpack.Target(_chain3(), count=2).with_name("peo")], max_loops=20)
        )
        assert grown.natoms == 6
        res = _box_gencan().run([molpack.Target.fixed_from(grown)], max_loops=5)
        assert res.converged is True
        assert res.natoms == 6
        assert np.array_equal(res.positions, grown.positions)
