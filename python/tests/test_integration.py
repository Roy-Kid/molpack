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


def _packer() -> molpack.Molpack:
    return molpack.Molpack().with_tolerance(2.0).with_progress(False)


@pytest.mark.integration
class TestSmallBoxPack:
    def test_single_species_converges(self, water_frame):
        target = (
            molpack.Target(water_frame, count=20)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
            )
        )
        result = _packer().with_seed(42).pack([target], max_loops=200)

        assert result.converged
        assert result.natoms == 60
        assert result.positions.shape == (60, 3)
        assert len(result.elements) == 60

    def test_two_species_pack(self, water_frame, urea_frame):
        box = molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [25.0, 25.0, 25.0])

        water = (
            molpack.Target(water_frame, count=20).with_name("water").with_restraint(box)
        )
        urea = (
            molpack.Target(urea_frame, count=10).with_name("urea").with_restraint(box)
        )

        result = _packer().with_seed(1_234_567).pack([water, urea], max_loops=200)

        w_natoms = water_frame["atoms"].nrows
        u_natoms = urea_frame["atoms"].nrows
        assert result.natoms == 20 * w_natoms + 10 * u_natoms
        assert result.positions.shape[0] == result.natoms


@pytest.mark.integration
class TestReproducibility:
    def test_same_seed_same_result(self, water_frame):
        def make() -> molpack.Target:
            return (
                molpack.Target(water_frame, count=10)
                .with_name("water")
                .with_restraint(
                    molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
                )
            )

        packer = molpack.Molpack().with_progress(False).with_seed(7)
        r1 = packer.pack([make()], max_loops=100)
        r2 = packer.pack([make()], max_loops=100)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_different_seeds_differ(self, water_frame):
        target = (
            molpack.Target(water_frame, count=10)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
            )
        )
        packer = molpack.Molpack().with_progress(False)
        r1 = packer.with_seed(1).pack([target], max_loops=100)
        r2 = packer.with_seed(2).pack([target], max_loops=100)
        assert not np.array_equal(r1.positions, r2.positions)


@pytest.mark.integration
class TestCompositeRestraints:
    def test_box_and_outside_sphere(self, water_frame):
        target = (
            molpack.Target(water_frame, count=5)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
            )
            .with_restraint(molpack.OutsideSphereRestraint([15.0, 15.0, 15.0], 5.0))
        )
        result = _packer().with_seed(42).pack([target], max_loops=150)

        centres = result.positions
        sphere_centre = np.array([15.0, 15.0, 15.0])
        dists = np.linalg.norm(centres - sphere_centre, axis=1)
        if result.converged:
            assert (dists >= 4.5).all()
