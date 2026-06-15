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


def _ion_frame() -> molrs.Frame:
    return molrs.Frame.from_dict(
        {
            "blocks": {
                "atoms": {
                    "x": np.array([0.0]),
                    "y": np.array([0.0]),
                    "z": np.array([0.0]),
                    "element": ["NA"],
                }
            }
        }
    )


def _packer() -> molpack.Molpack:
    return molpack.Molpack().with_progress(False)


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


class TestPackingBehavior:
    def test_native_plane_reproduces_gaussian(self):
        """GaussianPlane drives a species onto its target slab band."""
        mu, sigma = 20.0, 4.0
        target = (
            molpack.Target(_ion_frame(), count=200)
            .with_name("NA")
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(
                molpack.GaussianPlane([0.0, 0.0, 1.0], 0.0, 1000.0, mu, sigma)
            )
        )
        result = (
            _packer()
            .with_seed(1)
            .with_tolerance(2.0)
            .pack_with_report([target], max_loops=60)
        )
        z = np.asarray(result.positions)[:, 2]
        assert abs(z.mean() - mu) < 1.0, f"mean {z.mean():.2f} != {mu}"
        assert abs(z.std() - sigma) < 1.0, f"std {z.std():.2f} != {sigma}"

    def test_native_point_reproduces_shell(self):
        """GaussianPoint drives a species onto a spherical shell of radius mu."""
        mu, sigma = 14.0, 2.0
        center = [20.0, 20.0, 20.0]
        target = (
            molpack.Target(_ion_frame(), count=200)
            .with_name("NA")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [40.0, 40.0, 40.0])
            )
            .with_restraint(molpack.GaussianPoint(center, 1000.0, mu, sigma))
        )
        result = (
            _packer()
            .with_seed(1)
            .with_tolerance(2.0)
            .pack_with_report([target], max_loops=80)
        )
        pos = np.asarray(result.positions)
        r = np.linalg.norm(pos - np.asarray(center), axis=1)
        assert abs(r.mean() - mu) < 1.5, f"shell radius mean {r.mean():.2f} != {mu}"
        assert abs(r.std() - sigma) < 1.5, (
            f"shell thickness std {r.std():.2f} != {sigma}"
        )

    def test_native_exponential_plane_decays_from_wall(self):
        """ExponentialPlane drives a species into a decaying layer at the wall:
        mean depth ~ lambda, and the layer is densest near z=0 (most ions below
        the mean) rather than uniform."""
        lam = 8.0
        target = (
            molpack.Target(_ion_frame(), count=200)
            .with_name("NA")
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(molpack.ExponentialPlane([0.0, 0.0, 1.0], 0.0, 1000.0, lam))
        )
        result = (
            _packer()
            .with_seed(1)
            .with_tolerance(2.0)
            .pack_with_report([target], max_loops=60)
        )
        z = np.asarray(result.positions)[:, 2]
        # Exponential mean is lambda; densest at the wall → median < mean.
        assert abs(z.mean() - lam) < 3.0, f"mean depth {z.mean():.2f} != ~{lam}"
        assert np.median(z) < z.mean(), "exponential layer should be right-skewed"

    def test_native_tabulated_reproduces_arbitrary_profile(self):
        """TabulatedPlane packs a species onto an arbitrary target density
        (here a Gouy-Chapman-shaped counter-ion profile), matching its quantiles."""
        z = np.linspace(0.0, 50.0, 400)
        kappa, gamma = 1.0 / 4.0, np.tanh(2.9 / 4)
        cp = ((1 + gamma * np.exp(-kappa * z)) / (1 - gamma * np.exp(-kappa * z))) ** 2
        target = (
            molpack.Target(_ion_frame(), count=300)
            .with_name("NA")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [16.0, 16.0, 50.0])
            )
            .with_restraint(
                molpack.TabulatedPlane(
                    [0.0, 0.0, 1.0], 0.0, 1000.0, z.tolist(), cp.tolist()
                )
            )
        )
        result = (
            _packer()
            .with_seed(3)
            .with_tolerance(2.0)
            .pack_with_report([target], max_loops=60)
        )
        zz = np.asarray(result.positions)[:, 2]
        # Target quantiles from the supplied grid; packed sorted z should match.
        cdf = np.concatenate([[0.0], np.cumsum(0.5 * (cp[1:] + cp[:-1]) * np.diff(z))])
        cdf /= cdf[-1]
        n = len(zz)
        qt = np.interp((np.arange(n) + 0.5) / n, cdf, z)
        rms = np.sqrt(np.mean((np.sort(zz) - qt) ** 2))
        assert rms < 1.0, (
            f"packed profile deviates from prior: quantile-RMS {rms:.3f} A"
        )

    def test_duck_typed_collective_pulls_mean(self):
        """A duck-typed collective restraint shifts the group mean toward its
        target (well below the box centre at z=20)."""
        target = 12.0
        tgt = (
            molpack.Target(_ion_frame(), count=60)
            .with_name("NA")
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(MeanTether([0.0, 0.0, 1.0], target, 5000.0))
        )
        result = (
            _packer()
            .with_seed(1)
            .with_tolerance(2.0)
            .pack_with_report([tgt], max_loops=60)
        )
        z = np.asarray(result.positions)[:, 2]
        assert z.mean() < 17.0, f"mean {z.mean():.2f} not pulled toward {target}"


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
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(Recorder())
        )
        _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
            [target], max_loops=20
        )

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
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(Explodes())
        )
        with pytest.raises(ValueError, match="boom from collective"):
            _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
                [target], max_loops=20
            )

    def test_wrong_gradient_count_is_reraised(self):
        class WrongLen:
            def f(self, coords, s, s2):
                return 0.0

            def fg(self, coords, s, s2):
                return 0.0, [(0.0, 0.0, 0.0)]  # too few gradients

        target = (
            molpack.Target(_ion_frame(), count=60)
            .with_restraint(molpack.InsideBoxRestraint(BOX_LO, BOX_HI))
            .with_restraint(WrongLen())
        )
        with pytest.raises(TypeError, match="gradients for"):
            _packer().with_seed(1).with_tolerance(2.0).pack_with_report(
                [target], max_loops=20
            )
