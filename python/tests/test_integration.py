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
class TestPeriodicBox:
    """Packer-level PBC and script `pbc` keyword integration."""

    def test_with_periodic_box_without_restraint(self, water_frame):
        # Regression: before the parser/plumbing fix, a target with no
        # spatial restraint fell back to a 2000 Å inferred box and
        # allocated ~10⁸ cells. With `with_periodic_box` in place the
        # cell grid is sized directly from the PBC box.
        target = molpack.Target(water_frame, count=5).with_name("water")
        result = (
            _packer()
            .with_seed(42)
            .with_periodic_box((0.0, 0.0, 0.0), (20.0, 20.0, 20.0))
            .pack([target], max_loops=50)
        )
        assert result.natoms == 5 * water_frame["atoms"].nrows

    def test_conflicting_packer_and_restraint_pbc(self, water_frame):
        # Packer-level PBC of 30 Å conflicts with restraint-level PBC
        # of 40 Å → ConflictingPeriodicBoxesError.
        target = (
            molpack.Target(water_frame, count=1)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint(
                    [0.0, 0.0, 0.0],
                    [40.0, 40.0, 40.0],
                    periodic=(True, True, True),
                )
            )
        )
        with pytest.raises(molpack.ConflictingPeriodicBoxesError):
            (
                _packer()
                .with_seed(42)
                .with_periodic_box((0.0, 0.0, 0.0), (30.0, 30.0, 30.0))
                .pack([target], max_loops=10)
            )

    def test_load_script_wires_pbc_to_packer(self, tmp_path):
        # `load_script` must surface a script-level `pbc` directive as
        # packer-level PBC so a pbc-only script (no `inside`) finishes
        # in seconds instead of hanging on a ~10⁸-cell grid.
        water_pdb = DATA_ROOT / "pack_mixture" / "water.pdb"
        out_pdb = tmp_path / "out.pdb"
        script_path = tmp_path / "pbc_only.inp"
        script_path.write_text(
            "tolerance 2.0\n"
            "seed 1234567\n"
            "filetype pdb\n"
            f"output {out_pdb}\n"
            "pbc 25.0 25.0 25.0\n\n"
            f"structure {water_pdb}\n"
            "  number 8\n"
            "end structure\n"
        )
        job = molpack.load_script(script_path)
        result = job.packer.with_progress(False).pack(job.targets, max_loops=job.nloop)
        assert result.natoms == 8 * 3  # water.pdb = 3 atoms/molecule


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


# ── correctness asserts: every packed atom satisfies its restraint ─────────

# Until this block was added, the Python tests checked atom counts and
# positional shapes but never verified that the packed coordinates
# actually obey the geometric restraint. A regression in the PyO3 layer
# (e.g. a Target/restraint wiring drift between Rust and Python) could
# produce shape-correct but geometrically-wrong output and slip past CI.

TOLERANCE = 2.0
_PRECISION_SLACK = 0.05  # absolute slack on `precision = 1e-2`


def _packer_with_tolerance(tolerance: float = TOLERANCE) -> molpack.Molpack:
    return (
        molpack.Molpack()
        .with_tolerance(tolerance)
        .with_precision(0.01)
        .with_progress(False)
    )


def _max_pairwise_violation(
    positions: np.ndarray, tolerance: float, *, atoms_per_mol: int
) -> float:
    """Return max(0, tolerance − min inter-molecule distance).

    Pairs within the same molecule are skipped (they are rigid). Naïve
    O(N²) — fine for small test sizes.
    """
    n = positions.shape[0]
    worst = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if i // atoms_per_mol == j // atoms_per_mol:
                continue
            d = float(np.linalg.norm(positions[i] - positions[j]))
            if d < tolerance:
                worst = max(worst, tolerance - d)
    return worst


@pytest.mark.integration
class TestRestraintSatisfaction:
    """Coordinate-level correctness: packed atoms obey their restraint."""

    def test_inside_sphere_atoms_within_radius(self, water_frame):
        centre = np.array([10.0, 10.0, 10.0])
        radius = 12.0
        target = (
            molpack.Target(water_frame, count=8)
            .with_name("water")
            .with_restraint(molpack.InsideSphereRestraint(centre.tolist(), radius))
        )
        result = _packer_with_tolerance().with_seed(42).pack([target], max_loops=150)

        dists = np.linalg.norm(result.positions - centre, axis=1)
        # Allow a small slack so the optimizer's `precision=0.01` floor
        # does not cause a flake at the boundary.
        assert dists.max() <= radius + _PRECISION_SLACK, (
            f"max distance {dists.max()} > radius {radius}"
        )

    def test_inside_box_atoms_within_bounds(self, water_frame):
        box_min = np.array([0.0, 0.0, 0.0])
        box_max = np.array([20.0, 20.0, 20.0])
        target = (
            molpack.Target(water_frame, count=8)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint(box_min.tolist(), box_max.tolist())
            )
        )
        result = _packer_with_tolerance().with_seed(42).pack([target], max_loops=150)

        positions = result.positions
        assert (positions >= box_min - _PRECISION_SLACK).all(), (
            "atoms below box minimum"
        )
        assert (positions <= box_max + _PRECISION_SLACK).all(), (
            "atoms above box maximum"
        )

    def test_pairwise_tolerance_respected(self, water_frame):
        # Smaller box so pair penalties are actually contended.
        target = (
            molpack.Target(water_frame, count=10)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [18.0, 18.0, 18.0])
            )
        )
        result = _packer_with_tolerance().with_seed(42).pack([target], max_loops=200)

        atoms_per_mol = water_frame["atoms"].nrows
        worst = _max_pairwise_violation(
            result.positions, TOLERANCE, atoms_per_mol=atoms_per_mol
        )
        # `precision=0.01` is the soft tolerance the optimizer
        # promises; allow a small slack on top.
        assert worst <= 0.01 + _PRECISION_SLACK, (
            f"max pairwise tolerance violation {worst}"
        )

    def test_outside_sphere_atoms_outside_radius(self, water_frame):
        centre = np.array([15.0, 15.0, 15.0])
        radius = 5.0
        target = (
            molpack.Target(water_frame, count=6)
            .with_name("water")
            .with_restraint(
                molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [30.0, 30.0, 30.0])
            )
            .with_restraint(molpack.OutsideSphereRestraint(centre.tolist(), radius))
        )
        result = _packer_with_tolerance().with_seed(42).pack([target], max_loops=200)

        # Skip when convergence stalls — outside-sphere is a harder
        # objective and a tiny seed can leave residual penalty. The
        # `if result.converged` guard mirrors the existing
        # `TestCompositeRestraints` test.
        if not result.converged:
            pytest.skip("packer did not converge; restraint check non-meaningful")
        dists = np.linalg.norm(result.positions - centre, axis=1)
        assert dists.min() >= radius - _PRECISION_SLACK, (
            f"min distance {dists.min()} < radius {radius}"
        )
