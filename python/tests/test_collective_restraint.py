"""Tests for collective (group-level) restraints.

Covers the group path of the unified ``Target.with_restraint``: the built-in
compiled distribution-matchers (``molpack.GaussianPlane`` / ``GaussianPoint``)
and duck-typed Python restraints whose ``f``/``fg`` receive every copy of a
species at once and return per-copy gradients.
"""

from __future__ import annotations

import molrs
import numpy as np
import pytest

import molpack

BOX_LO = [0.0, 0.0, 0.0]
BOX_HI = [12.0, 12.0, 40.0]


def _ion_frame() -> molrs.store.Frame:
    return molrs.store.Frame(
        {
            "atoms": {
                "x": np.array([0.0]),
                "y": np.array([0.0]),
                "z": np.array([0.0]),
                "element": ["NA"],
            }
        }
    )


def _packer() -> molpack.GenCanPack:
    return molpack.GenCanPack().with_progress(False)


class MeanTether:
    """Minimal collective restraint: pull the group's *mean* coordinate along a
    plane normal toward ``target``. ``L = (lam/2)(mean(xi) - target)^2`` couples
    every copy through the shared mean, so it exercises the group-level contract
    without any external dependency."""

    def __init__(self, normal, target, strength):
        n = np.asarray(normal, dtype=np.float64)
        self.n = n / np.linalg.norm(n)
        self.target = float(target)
        self.lam = float(strength)

    def _xi(self, coords):
        return np.asarray(coords, dtype=np.float64) @ self.n

    def f(self, coords, scale, scale2):
        d = self._xi(coords).mean() - self.target
        return 0.5 * self.lam * d * d

    def fg(self, coords, scale, scale2):
        xi = self._xi(coords)
        n = len(xi)
        d = xi.mean() - self.target
        g = ((self.lam * d / n) * np.ones(n))[:, None] * self.n[None, :]
        return 0.5 * self.lam * d * d, g.tolist()


class TestAttachment:
    def test_native_plane_accepted(self):
        t = molpack.Target(_ion_frame(), count=4).with_restraint(
            molpack.GaussianPlane([0.0, 0.0, 1.0], 0.0, 100.0, 20.0, 4.0)
        )
        assert t is not None

    def test_native_point_accepted(self):
        t = molpack.Target(_ion_frame(), count=4).with_restraint(
            molpack.GaussianPoint([0.0, 0.0, 0.0], 100.0, 20.0, 4.0)
        )
        assert t is not None

    def test_native_exponential_accepted(self):
        for r in (
            molpack.ExponentialPlane([0.0, 0.0, 1.0], 0.0, 100.0, 5.0),
            molpack.ExponentialPoint([0.0, 0.0, 0.0], 100.0, 5.0),
        ):
            t = molpack.Target(_ion_frame(), count=4).with_restraint(r)
            assert t is not None

    def test_exponential_rejects_nonpositive_lambda(self):
        with pytest.raises(ValueError, match="lambda"):
            molpack.ExponentialPlane([0.0, 0.0, 1.0], 0.0, 100.0, 0.0)
        with pytest.raises(ValueError, match="lambda"):
            molpack.ExponentialPoint([0.0, 0.0, 0.0], 100.0, -1.0)

    def test_native_tabulated_accepted(self):
        xs = [0.0, 1.0, 2.0, 3.0]
        rho = [4.0, 2.0, 1.0, 0.5]
        for r in (
            molpack.TabulatedPlane([0.0, 0.0, 1.0], 0.0, 100.0, xs, rho),
            molpack.TabulatedPoint([0.0, 0.0, 0.0], 100.0, xs, rho),
        ):
            t = molpack.Target(_ion_frame(), count=4).with_restraint(r)
            assert t is not None

    def test_tabulated_rejects_bad_grid(self):
        with pytest.raises(ValueError, match="ascending"):
            molpack.TabulatedPlane([0.0, 0.0, 1.0], 0.0, 100.0, [1.0, 0.0], [1.0, 1.0])
        with pytest.raises(ValueError, match="positive total mass"):
            molpack.TabulatedPlane([0.0, 0.0, 1.0], 0.0, 100.0, [0.0, 1.0], [0.0, 0.0])

    def test_duck_typed_accepted(self):
        t = molpack.Target(_ion_frame(), count=4).with_restraint(
            MeanTether([0.0, 0.0, 1.0], 20.0, 100.0)
        )
        assert t is not None

    def test_object_without_methods_rejected(self):
        class Empty:
            pass

        with pytest.raises(TypeError, match="expected a restraint"):
            molpack.Target(_ion_frame(), count=1).with_restraint(Empty())

    def test_plane_rejects_nonpositive_sigma(self):
        with pytest.raises(ValueError, match="sigma"):
            molpack.GaussianPlane([0.0, 0.0, 1.0], 0.0, 100.0, 20.0, 0.0)

    def test_plane_rejects_zero_normal(self):
        with pytest.raises(ValueError, match="normal"):
            molpack.GaussianPlane([0.0, 0.0, 0.0], 0.0, 100.0, 20.0, 4.0)

    def test_point_rejects_nonpositive_sigma(self):
        with pytest.raises(ValueError, match="sigma"):
            molpack.GaussianPoint([0.0, 0.0, 0.0], 100.0, 20.0, 0.0)

    def test_self_separation_accepted(self):
        t = molpack.Target(_ion_frame(), count=4).with_restraint(
            molpack.SelfSeparation(10.0)
        )
        assert t is not None

    def test_self_separation_strength_defaults_to_one(self):
        r = molpack.SelfSeparation(10.0)
        assert r.d_min == 10.0
        assert "10" in repr(r)

    def test_self_separation_rejects_nonpositive_distance(self):
        with pytest.raises(ValueError, match="d_min"):
            molpack.SelfSeparation(0.0)
        with pytest.raises(ValueError, match="d_min"):
            molpack.SelfSeparation(-1.0)

    def test_self_separation_rejects_nonpositive_strength(self):
        with pytest.raises(ValueError, match="strength"):
            molpack.SelfSeparation(10.0, 0.0)


class TestCallContract:
    def test_fg_receives_whole_group(self):
        # A non-trivial copy count forces the main objective phase, where the
        # collective term runs (it is gated off during independent placement).
        count = 60
        seen: dict = {}

        class Recorder:
            def f(self, coords, scale, scale2):
                return 0.0

            def fg(self, coords, scale, scale2):
                seen["coords"] = coords
                seen["scale"] = scale
                return 0.0, [(0.0, 0.0, 0.0)] * len(coords)

        target = (
            molpack.Target(_ion_frame(), count=count)
            .with_restraint(molrs.spatial.Cuboid(BOX_LO, np.subtract(BOX_HI, BOX_LO)))
            .with_restraint(Recorder())
        )
        _packer().with_seed(1).with_tolerance(2.0).run([target], max_loops=20)

        assert "coords" in seen, "collective fg was never called"
        coords = seen["coords"]
        assert len(coords) == count, "fg must receive every copy at once"
        assert len(coords[0]) == 3
        assert isinstance(seen["scale"], float)


class TestErrorPropagation:
    def test_exception_in_fg_is_reraised(self):
        class Explodes:
            def f(self, coords, s, s2):
                return 0.0

            def fg(self, coords, s, s2):
                raise ValueError("boom from collective")

        target = (
            molpack.Target(_ion_frame(), count=60)
            .with_restraint(molrs.spatial.Cuboid(BOX_LO, np.subtract(BOX_HI, BOX_LO)))
            .with_restraint(Explodes())
        )
        with pytest.raises(ValueError, match="boom from collective"):
            _packer().with_seed(1).with_tolerance(2.0).run([target], max_loops=20)

    def test_wrong_gradient_count_is_reraised(self):
        class WrongLen:
            def f(self, coords, s, s2):
                return 0.0

            def fg(self, coords, s, s2):
                return 0.0, [(0.0, 0.0, 0.0)]  # too few gradients

        target = (
            molpack.Target(_ion_frame(), count=60)
            .with_restraint(molrs.spatial.Cuboid(BOX_LO, np.subtract(BOX_HI, BOX_LO)))
            .with_restraint(WrongLen())
        )
        with pytest.raises(TypeError, match="gradients for"):
            _packer().with_seed(1).with_tolerance(2.0).run([target], max_loops=20)


class TestSelfSeparation:
    """The anti-clustering restraint: copies of one species keep their distance
    from each other, which nothing in the pair term asks of them."""

    BOX = ([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])  # origin, lengths
    N = 27
    D_MIN = 10.0

    def _run(self, *, separate: bool, seed: int = 1):
        target = (
            molpack.Target(_ion_frame(), count=self.N)
            .with_name("NA")
            .with_restraint(molrs.spatial.Cuboid(*self.BOX))
        )
        if separate:
            target = target.with_restraint(molpack.SelfSeparation(self.D_MIN))
        result = (
            _packer().with_seed(seed).with_tolerance(2.0).run([target], max_loops=200)
        )
        return np.asarray(result.positions), result

    @staticmethod
    def _min_pair_distance(pos):
        d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        return d.min()

    def test_an_impossible_request_is_reported_not_hidden(self):
        # 27 copies cannot be 30 Å apart in a 40 Å box. The run must not claim
        # success: the penalty stays in frest and convergence fails.
        target = (
            molpack.Target(_ion_frame(), count=self.N)
            .with_name("NA")
            .with_restraint(molrs.spatial.Cuboid(*self.BOX))
            .with_restraint(molpack.SelfSeparation(30.0))
        )
        result = _packer().with_seed(1).with_tolerance(2.0).run([target], max_loops=40)
        assert not result.converged
        assert result.frest > 0.0

    def test_duck_typed_restraints_still_receive_scale_arguments(self):
        # The Rust seam grew a context argument; the Python contract must not
        # have. A duck-typed collective restraint still gets (coords, scale,
        # scale2) and must compose with a native one on the same species.
        target = (
            molpack.Target(_ion_frame(), count=20)
            .with_name("NA")
            .with_restraint(molrs.spatial.Cuboid(BOX_LO, np.subtract(BOX_HI, BOX_LO)))
            .with_restraint(MeanTether([0.0, 0.0, 1.0], 20.0, 100.0))
            .with_restraint(molpack.SelfSeparation(4.0))
        )
        result = _packer().with_seed(1).with_tolerance(2.0).run([target], max_loops=60)
        z = np.asarray(result.positions)[:, 2]
        assert abs(z.mean() - 20.0) < 3.0
