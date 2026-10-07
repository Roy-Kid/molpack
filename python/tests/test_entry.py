"""Engine entries (engine-entry-split): GenCanPack / CbmcGrow 1:1 bindings."""

from __future__ import annotations

import molrs
import numpy as np
import pytest

import molpack


def _dimer():

    fr = molrs.core.Frame()
    fr["atoms"] = {
        "x": np.array([0.0, 1.5]),
        "y": np.zeros(2),
        "z": np.zeros(2),
        "element": ["C", "C"],
    }
    return fr


def _chain5():

    n = 5
    fr = molrs.core.Frame()
    fr["atoms"] = {
        "x": np.arange(n) * 1.5,
        "y": np.zeros(n),
        "z": np.zeros(n),
        "element": ["C"] * n,
    }
    fr["bonds"] = {"atomi": np.arange(n - 1), "atomj": np.arange(1, n)}
    return fr


class TestGenCanPackEntry:
    def test_public_surface_is_state_not_pack_result(self):
        assert hasattr(molpack, "State")
        assert not hasattr(molpack, "PackResult")
        engine = molpack.GenCanPack()
        assert hasattr(engine, "with_restart")
        assert not hasattr(engine, "seeded_from")

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
        assert res.degraded == 0

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
            .with_restart(grown)
            .with_seed(9)
            .with_tolerance(1.0)
            .run([target()], max_loops=60)
        )
        assert pushed.natoms == grown.natoms
        assert pushed.degraded == 0
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
                .with_restart(grown)
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


class TestEngineErrorPaths:
    def test_empty_targets_list_raises(self):
        packer = molpack.GenCanPack().with_progress(False).with_seed(1)
        with pytest.raises(molpack.NoTargetsError):
            packer.run([], max_loops=10)

    def test_invalid_pbc_raises_typed_error(self):

        positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        frame = molrs.core.Frame(
            {
                "atoms": {
                    "x": positions[:, 0],
                    "y": positions[:, 1],
                    "z": positions[:, 2],
                    "element": ["X"],
                }
            }
        )
        target = molpack.Target(frame, 1)
        packer = (
            molpack.GenCanPack()
            .with_progress(False)
            .with_seed(1)
            .with_periodic_box((0.0, 0.0, 0.0), (0.0, 10.0, 10.0))
        )
        with pytest.raises(molpack.InvalidPBCBoxError):
            packer.run([target], max_loops=10)

    def test_pack_error_is_runtime_error_subclass(self):
        assert issubclass(molpack.NoTargetsError, molpack.PackError)
        assert issubclass(molpack.PackError, RuntimeError)
