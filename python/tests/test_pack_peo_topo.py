"""Topological PEO templates from molrs/molpy, packed by molpack.

Mirrors ``python/examples/pack_peo_topo.py``: SMILES + conformer +
``PolymerBuilder``, then the packing named-error and tiny-run contracts.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

import molpack

_EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "pack_peo_topo.py"
_spec = importlib.util.spec_from_file_location("pack_peo_topo", _EXAMPLE)
assert _spec is not None and _spec.loader is not None
topo = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(topo)


class TestTemplates:
    def test_star_is_a_tree(self):
        star = topo.make_star(2, seed=42)
        n_at = star.n_atoms
        n_bd = len(list(star.bonds))
        assert n_bd == n_at - 1
        assert n_at == 77  # 4-arm, arm_length=2, X4 core; smoke 2026-09-04

    def test_ring_is_unicyclic(self):
        ring = topo.make_ring(4, seed=42)
        n_at = ring.n_atoms
        n_bd = len(list(ring.bonds))
        assert n_bd == n_at
        assert n_at == 28

    def test_star_arm_length_rejects_zero(self):
        with pytest.raises(ValueError, match="arm_length"):
            topo.make_star(0)

    def test_ring_rejects_too_small(self):
        with pytest.raises(ValueError, match="n >= 3"):
            topo.make_ring(2)

    def test_linear_is_a_tree(self):
        linear = topo.make_linear(2, seed=42)
        n_at = linear.n_atoms
        n_bd = len(list(linear.bonds))
        assert n_bd == n_at - 1
        assert n_at == 17  # two EO residues after condensation; smoke 2026-09-05

    def test_linear_n_rejects_zero(self):
        with pytest.raises(ValueError, match="linear n"):
            topo.make_linear(0)


class TestNamedErrors:
    def test_cbmc_refuses_ring(self):
        ring = topo.make_ring(4, seed=1)
        target = molpack.Target(ring.to_frame(), 1).with_name("c-PEO")
        prior = molpack.TorsionPrior.three_state_from_c_inf(topo.PEO_C_INF, topo.TET)
        with pytest.raises(ValueError, match="ring"):
            molpack.CbmcGrow(prior).with_seed(1).with_tolerance(2.0).with_density(
                0.2
            ).with_progress(False).run([target], max_loops=4)

    def test_lattice_refuses_ring(self):
        ring = topo.make_ring(4, seed=1)
        target = molpack.Target(ring.to_frame(), 1).with_name("c-PEO")
        prior = molpack.TorsionPrior.three_state_from_c_inf(topo.PEO_C_INF, topo.TET)
        with pytest.raises(ValueError, match="ring"):
            molpack.LatticeGrow(prior).with_seed(1).with_tolerance(2.0).with_density(
                0.2
            ).with_progress(False).run([target], max_loops=4)


class TestTinyPack:
    def test_lattice_grows_star(self):
        star = topo.make_star(2, seed=42)
        target = molpack.Target(star.to_frame(), 1).with_name("star-PEO")
        prior = molpack.TorsionPrior.three_state_from_c_inf(topo.PEO_C_INF, topo.TET)
        grown = (
            molpack.LatticeGrow(prior)
            .with_seed(42)
            .with_tolerance(2.0)
            .with_density(0.2)
            .with_progress(False)
            .run([target], max_loops=40)
        )
        assert grown.natoms == star.n_atoms

    def test_cbmc_grows_one_star(self):
        star = topo.make_star(2, seed=42)
        target = molpack.Target(star.to_frame(), 1).with_name("star-PEO")
        prior = molpack.TorsionPrior.three_state_from_c_inf(topo.PEO_C_INF, topo.TET)
        grown = (
            molpack.CbmcGrow(prior)
            .with_seed(42)
            .with_tolerance(2.0)
            .with_density(0.2)
            .with_progress(False)
            .run([target], max_loops=40)
        )
        assert grown.natoms == star.n_atoms
        assert grown.converged
        assert grown.degraded == 0
        assert grown.fdist == 0.0

    def test_auhl_two_stars_push_off_converges(self):
        star = topo.make_star(2, seed=42)
        target = molpack.Target(star.to_frame(), 2).with_name("star-PEO")
        prior = molpack.TorsionPrior.three_state_from_c_inf(topo.PEO_C_INF, topo.TET)
        grown = (
            molpack.CbmcGrow(prior)
            .with_seed(42)
            .with_tolerance(0.6)
            .with_density(0.5)
            .with_progress(False)
            .run([target], max_loops=40)
        )
        pushed = (
            molpack.GenCanPack()
            .with_restart(grown)
            .with_seed(42)
            .with_tolerance(2.0)
            .with_progress(False)
            .run([target], max_loops=80)
        )
        assert pushed.natoms == 2 * star.n_atoms
        assert pushed.converged
        assert pushed.fdist == 0.0

    def test_gencan_packs_two_rings(self):
        ring = topo.make_ring(4, seed=42)
        target = molpack.Target(ring.to_frame(), 2).with_name("c-PEO")
        packed = (
            molpack.GenCanPack()
            .with_seed(42)
            .with_tolerance(2.0)
            .with_density(0.3)
            .with_progress(False)
            .run([target], max_loops=40)
        )
        assert packed.natoms == 2 * ring.n_atoms
        assert packed.converged
        assert packed.fdist == 0.0
