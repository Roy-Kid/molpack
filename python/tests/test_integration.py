"""End-to-end integration tests using real PDB fixtures.

Exercises the full loader → Target → Molpack → PackResult pipeline at
small scale. For the full Packmol-equivalent workloads see the
scripts under ``python/examples/``.
"""

from __future__ import annotations

from pathlib import Path

import molrs
import numpy as np
import pytest

import molpack

DATA_ROOT = Path(__file__).resolve().parent.parent.parent / "examples"


@pytest.fixture(scope="module")
def water_frame():
    return molrs.read_pdb(str(DATA_ROOT / "pack_mixture" / "water.pdb"))


@pytest.fixture(scope="module")
def urea_frame():
    return molrs.read_pdb(str(DATA_ROOT / "pack_mixture" / "urea.pdb"))


@pytest.mark.integration
class TestSmallBoxPack:
    def test_single_species_converges(self, water_frame):
        target = molpack.Target("water", water_frame, count=20).with_restraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        )
        packer = molpack.Molpack(tolerance=2.0, progress=False)
        result = packer.pack([target], max_loops=200, seed=42)

        assert result.converged
        assert result.natoms == 60
        assert result.positions.shape == (60, 3)
        assert len(result.elements) == 60

    def test_two_species_pack(self, water_frame, urea_frame):
        box = molpack.InsideBox([0.0, 0.0, 0.0], [25.0, 25.0, 25.0])

        water = molpack.Target("water", water_frame, count=20).with_restraint(box)
        urea = molpack.Target("urea", urea_frame, count=10).with_restraint(box)

        packer = molpack.Molpack(tolerance=2.0, progress=False)
        result = packer.pack([water, urea], max_loops=200, seed=1_234_567)

        w_natoms = water_frame["atoms"].nrows
        u_natoms = urea_frame["atoms"].nrows
        assert result.natoms == 20 * w_natoms + 10 * u_natoms
        assert result.positions.shape[0] == result.natoms


@pytest.mark.integration
class TestReproducibility:
    def test_same_seed_same_result(self, water_frame):
        make = lambda: (  # noqa: E731
            molpack.Target("water", water_frame, count=10).with_restraint(
                molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
            )
        )
        packer = molpack.Molpack(progress=False)
        r1 = packer.pack([make()], max_loops=100, seed=7)
        r2 = packer.pack([make()], max_loops=100, seed=7)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_different_seeds_differ(self, water_frame):
        target = molpack.Target("water", water_frame, count=10).with_restraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
        )
        packer = molpack.Molpack(progress=False)
        r1 = packer.pack([target], max_loops=100, seed=1)
        r2 = packer.pack([target], max_loops=100, seed=2)
        assert not np.array_equal(r1.positions, r2.positions)


@pytest.mark.integration
class TestCompositeRestraints:
    def test_box_and_outside_sphere(self, water_frame):
        target = (
            molpack.Target("water", water_frame, count=5)
            .with_restraint(molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0]))
            .with_restraint(molpack.OutsideSphere(5.0, [15.0, 15.0, 15.0]))
        )
        packer = molpack.Molpack(tolerance=2.0, progress=False)
        result = packer.pack([target], max_loops=150, seed=42)

        centres = result.positions
        sphere_centre = np.array([15.0, 15.0, 15.0])
        dists = np.linalg.norm(centres - sphere_centre, axis=1)
        if result.converged:
            assert (dists >= 4.5).all()
