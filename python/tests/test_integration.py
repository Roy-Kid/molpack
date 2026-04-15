"""End-to-end integration tests using real PDB fixtures.

Exercises the full loader → Target → Packer → PackResult pipeline at
small scale. For the full Packmol-equivalent workloads see the
scripts under ``python/examples/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

import molpack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "examples"))

from _common import read_pdb_as_arrays  # noqa: E402

DATA_ROOT = HERE.parent.parent / "examples"


@pytest.fixture(scope="module")
def water_template():
    return read_pdb_as_arrays(DATA_ROOT / "pack_mixture" / "water.pdb")


@pytest.fixture(scope="module")
def urea_template():
    return read_pdb_as_arrays(DATA_ROOT / "pack_mixture" / "urea.pdb")


@pytest.mark.integration
class TestSmallBoxPack:
    def test_single_species_converges(self, water_template):
        pos, rad, els = water_template
        target = (
            molpack.Target.from_coords(pos, rad, count=20, elements=els)
            .with_name("water")
            .with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0]))
        )
        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=200, seed=42)

        assert result.converged
        assert result.natoms == 60
        assert result.positions.shape == (60, 3)
        assert len(result.elements) == 60

    def test_two_species_pack(self, water_template, urea_template):
        w_pos, w_rad, w_els = water_template
        u_pos, u_rad, u_els = urea_template

        box = molpack.InsideBox([0.0, 0.0, 0.0], [25.0, 25.0, 25.0])

        water = (
            molpack.Target.from_coords(w_pos, w_rad, count=20, elements=w_els)
            .with_name("water")
            .with_constraint(box)
        )
        urea = (
            molpack.Target.from_coords(u_pos, u_rad, count=10, elements=u_els)
            .with_name("urea")
            .with_constraint(box)
        )

        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([water, urea], max_loops=200, seed=1_234_567)

        assert result.natoms == 20 * 3 + 10 * 8
        assert result.positions.shape[0] == result.natoms
        # Packed elements should contain every input element.
        assert set(result.elements) == set(w_els) | set(u_els)


@pytest.mark.integration
class TestReproducibility:
    def test_same_seed_same_result(self, water_template):
        pos, rad, els = water_template
        make = lambda: (  # noqa: E731
            molpack.Target.from_coords(
                pos, rad, count=10, elements=els
            ).with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0]))
        )
        packer = molpack.Packer().with_progress(False)
        r1 = packer.pack([make()], max_loops=100, seed=7)
        r2 = packer.pack([make()], max_loops=100, seed=7)
        np.testing.assert_array_equal(r1.positions, r2.positions)

    def test_different_seeds_differ(self, water_template):
        pos, rad, els = water_template
        target = molpack.Target.from_coords(
            pos, rad, count=10, elements=els
        ).with_constraint(molpack.InsideBox([0.0, 0.0, 0.0], [15.0, 15.0, 15.0]))
        packer = molpack.Packer().with_progress(False)
        r1 = packer.pack([target], max_loops=100, seed=1)
        r2 = packer.pack([target], max_loops=100, seed=2)
        assert not np.array_equal(r1.positions, r2.positions)


@pytest.mark.integration
class TestCompositeConstraints:
    def test_box_and_outside_sphere(self, water_template):
        pos, rad, els = water_template
        target = molpack.Target.from_coords(
            pos, rad, count=5, elements=els
        ).with_constraint(
            molpack.InsideBox([0.0, 0.0, 0.0], [30.0, 30.0, 30.0]).and_(
                molpack.OutsideSphere(5.0, [15.0, 15.0, 15.0])
            )
        )
        packer = molpack.Packer(tolerance=2.0).with_progress(False)
        result = packer.pack([target], max_loops=150, seed=42)

        # After packing the shell region, no atom should be inside the sphere.
        centres = result.positions
        sphere_centre = np.array([15.0, 15.0, 15.0])
        dists = np.linalg.norm(centres - sphere_centre, axis=1)
        # Allow for tolerance: converged run should place atoms outside r=5.
        if result.converged:
            assert (dists >= 4.5).all()
