"""Engine entries (engine-entry-split): GenCanPack / CbmcGrow 1:1 bindings."""

from __future__ import annotations

import numpy as np
import pytest

import molpack


def _dimer():
    import molrs

    fr = molrs.Frame()
    fr["atoms"] = {
        "x": np.array([0.0, 1.5]),
        "y": np.zeros(2),
        "z": np.zeros(2),
        "element": ["C", "C"],
    }
    return fr


def _chain5():
    import molrs

    n = 5
    fr = molrs.Frame()
    fr["atoms"] = {
        "x": np.arange(n) * 1.5,
        "y": np.zeros(n),
        "z": np.zeros(n),
        "element": ["C"] * n,
    }
    fr["bonds"] = {"atomi": np.arange(n - 1), "atomj": np.arange(1, n)}
    return fr


class TestGenCanPackEntry:
    def test_deterministic_under_seed(self):
        # Same seed, same targets → bitwise-identical positions. (Parity
        # against the deleted legacy `Molpack` was proven before its
        # removal — engine-entry-split migration record.)
        def pack():
            return (
                molpack.GenCanPack()
                .with_seed(11)
                .with_tolerance(2.0)
                .with_periodic_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
                .run([molpack.Target(_dimer(), count=6)], max_loops=50)
            )

        a = pack()
        b = pack()
        assert a.converged
        assert np.array_equal(a.positions, b.positions)

    def test_one_engine_one_run(self):
        eng = molpack.GenCanPack().with_periodic_box([0.0] * 3, [20.0] * 3)
        eng.run([molpack.Target(_dimer(), count=2)], max_loops=10)
        with pytest.raises(RuntimeError, match="one engine, one run"):
            eng.run([molpack.Target(_dimer(), count=2)], max_loops=10)


class TestCbmcGrowEntry:
    def test_grows_and_reports_honestly(self):
        res = (
            molpack.CbmcGrow(molpack.TorsionPrior.uniform())
            .with_seed(9)
            .with_tolerance(1.0)
            .with_periodic_box([0.0] * 3, [20.0] * 3)
            .run([molpack.Target(_chain5(), count=2)], max_loops=60)
        )
        assert res.natoms == 10
        assert res.converged is True
        assert res.softened == 0

    def test_seeded_push_off_chain(self):
        # The explicit push-off chain: grow, then continue the SAME free
        # targets on the grown state with a seeded GenCanPack.
        target = lambda: molpack.Target(_chain5(), count=2)  # noqa: E731
        grown = (
            molpack.CbmcGrow(molpack.TorsionPrior.uniform())
            .with_seed(9)
            .with_tolerance(1.0)
            .with_periodic_box([0.0] * 3, [20.0] * 3)
            .run([target()], max_loops=60)
        )
        pushed = (
            molpack.GenCanPack()
            .seeded_from(grown)
            .with_seed(9)
            .with_tolerance(1.0)
            .run([target()], max_loops=60)
        )
        assert pushed.natoms == grown.natoms
        assert pushed.softened == 0
        assert pushed.converged is True

    def test_seeded_shape_mismatch_is_named(self):
        grown = (
            molpack.CbmcGrow(molpack.TorsionPrior.uniform())
            .with_seed(9)
            .with_tolerance(1.0)
            .with_periodic_box([0.0] * 3, [20.0] * 3)
            .run([molpack.Target(_chain5(), count=2)], max_loops=60)
        )
        with pytest.raises(ValueError, match="seeded run"):
            (
                molpack.GenCanPack()
                .seeded_from(grown)
                .run([molpack.Target(_chain5(), count=3)], max_loops=10)
            )

    def test_chaining_over_fixed_matrix(self):
        grown = (
            molpack.CbmcGrow(molpack.TorsionPrior.uniform())
            .with_seed(9)
            .with_tolerance(1.0)
            .with_periodic_box([0.0] * 3, [20.0] * 3)
            .run([molpack.Target(_chain5(), count=2)], max_loops=60)
        )
        packed = (
            molpack.GenCanPack()
            .with_seed(3)
            .with_tolerance(1.0)
            .with_periodic_box([0.0] * 3, [20.0] * 3)
            .run(
                [molpack.Target.fixed_from(grown), molpack.Target(_dimer(), count=4)],
                max_loops=80,
            )
        )
        assert packed.converged is True
        assert packed.natoms == grown.natoms + 8
        assert np.array_equal(packed.positions[: grown.natoms], grown.positions)
