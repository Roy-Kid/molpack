"""Chain-growth entry bindings (``CbmcGrow``, engine-entry-split).

Covers the typed prior surface (``TorsionPrior`` / ``AnglePrior``), the
growth knobs on the entry itself, the named-error contracts (an unsupported
combination is a ``ValueError``, never a silent fall-back), the
``with_density`` mutual exclusion, and end-to-end grow packs from a real
bonded ``molrs`` frame.

Mirrors ``tests/grow.rs`` at the binding level: frames are genuine molrs
frames (never dicts or mocks), and bonded chains use the same planar-zigzag
template as the Rust suite's ``chain_frame`` fixture.
"""

from __future__ import annotations

import math

import molrs
import numpy as np
import pytest

from molpack import (
    AnglePrior,
    CbmcGrow,
    GenCanPack,
    LatticeGrow,
    Target,
    TorsionPrior,
)

#: CODATA Avogadro constant, exactly as fixed by the 2019 SI (the constant
#: inside ``with_density``'s box formula).
AVOGADRO = 6.02214076e23


def _chain_frame(n: int, bond: float = 1.53, bonds: bool = True) -> molrs.Frame:
    """Planar zigzag bead chain with tetrahedral (109.5°) angles — the Rust
    suite's ``chain_frame`` fixture. ``bonds=False`` drops the bonds block:
    the "bare coordinates" shape a grow target must reject by name."""
    theta = math.radians(109.5)
    alpha = (math.pi - theta) / 2.0
    dx, dz = bond * math.cos(alpha), bond * math.sin(alpha)
    idx = np.arange(n)
    blocks: dict = {
        "atoms": {
            "x": idx.astype(np.float64) * dx,
            "y": np.zeros(n, dtype=np.float64),
            "z": np.where(idx % 2 == 1, dz, 0.0),
            "element": ["C"] * n,
        }
    }
    if bonds:
        blocks["bonds"] = {
            "atomi": np.arange(0, n - 1, dtype=np.uint64),
            "atomj": np.arange(1, n, dtype=np.uint64),
        }
    return molrs.Frame(blocks)


def _grow() -> CbmcGrow:
    return CbmcGrow(TorsionPrior.uniform()).with_progress(False)


class TestTypedSurface:
    """The typed prior objects and entry builders — never strings."""

    def test_torsion_prior_constructors(self):
        # All four constructors build, and the repr names the variant — the
        # only readable window on a write-only (frozen) prior.
        assert "Uniform" in repr(TorsionPrior.uniform())
        assert "Template" in repr(TorsionPrior.template(2.0))
        assert "States" in repr(
            TorsionPrior.states([(0.0, 0.645), (2.094, 0.1775), (-2.094, 0.1775)])
        )
        # The C∞ calibration returns a discrete-states (RIS) prior (spec §5.1).
        calibrated = TorsionPrior.three_state_from_c_inf(5.5, math.radians(109.47))
        assert "States" in repr(calibrated)

    def test_angle_prior_constructors(self):
        assert "Template" in repr(AnglePrior.template())
        assert "Wlc" in repr(AnglePrior.wlc(3.0))
        assert "Wlc" in repr(AnglePrior.wlc_from_c_inf(1.76))

    def test_grow_entry_requires_torsion_prior(self):
        # The torsion prior is mandatory — no default, no empty constructor
        # (uniform sampling is quantitatively wrong for melts, spec §5.1).
        with pytest.raises(TypeError):
            CbmcGrow()  # ty: ignore[missing-argument]

    def test_grow_entry_builder_chain(self):
        base = CbmcGrow(TorsionPrior.uniform())
        chained = (
            base.with_trials(16)
            .with_selectivity(1.5)
            .with_soft_shell(0.8)
            .with_retract(6)
            .with_relax(20, 4)
            .with_soften_after(30)
            .with_min_hard_scale(0.85)
            .with_angle_prior(AnglePrior.template())
        )
        assert isinstance(chained, CbmcGrow)
        # Builders return a NEW entry; the original stays buildable.
        assert chained is not base
        assert repr(chained)

    def test_grow_has_no_exclusion_depth_knob(self):
        # Spec 05 unhooked the engine knob; the table lives on Target.
        entry = CbmcGrow(TorsionPrior.uniform())
        assert hasattr(entry, "with_exclusion_depth") is False


class TestNamedErrors:
    """Unsupported combinations are named ``ValueError``s — principle 3."""

    def test_grow_without_bonds_cannot_be_grown(self):
        # Bare coordinates (no bonds block) carry no chemistry to grow from;
        # molpack refuses by name instead of silently packing rigid bodies.
        target = Target(_chain_frame(5, bonds=False), 2)
        engine = (
            _grow()
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )
        with pytest.raises(ValueError, match="cannot be grown"):
            engine.run([target], max_loops=10)

    def test_grow_rejects_non_binary_special_bond(self):
        # Fractional 1-4 stores on Target; CbmcGrow.run wraps
        # GrowError::NonBinarySpecialBond in PackError::Grow (ValueError).
        target = Target(_chain_frame(5), 2).with_special_bonds([0.0, 0.0, 0.5, 1.0])
        engine = (
            _grow()
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )
        with pytest.raises(ValueError, match="1-4") as excinfo:
            engine.run([target], max_loops=10)
        assert "with_atom_radius" in str(excinfo.value)

    def test_grow_needs_box(self):
        # Growth needs the final volume from atom 0: no box, no cell, no
        # density → named error, not a guessed default.
        target = Target(_chain_frame(5), 2)
        engine = _grow().with_seed(7).with_tolerance(2.0)
        with pytest.raises(ValueError, match="needs a box"):
            engine.run([target], max_loops=10)

    def test_density_conflicts_with_box(self):
        # One source of truth for the volume — never a silent precedence rule.
        target = Target(_chain_frame(5), 2).with_mass(72.0)
        engine = (
            _grow()
            .with_seed(3)
            .with_tolerance(2.0)
            .with_density(0.05)
            .with_periodic_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        )
        with pytest.raises(ValueError, match="mutually exclusive"):
            engine.run([target], max_loops=10)


class TestGrowPack:
    """End-to-end growth through the wheel: same constructive guarantees as
    the Rust suite (``fdist == 0.0`` strict — hard rejection, not descent)."""

    def test_grow_pack_constructive(self):
        n, copies, edge = 6, 4, 20.0
        target = Target(_chain_frame(n), copies)
        result = (
            _grow()
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [edge, edge, edge])
            .run([target], max_loops=50)
        )
        assert result.converged
        assert result.softened == 0
        # Strict zero: hard-core violation is rejected during growth, never
        # penalized afterwards (ac-004).
        assert result.fdist == 0.0
        assert result.frest == 0.0

        pos = result.positions
        assert pos.shape == (copies * n, 3)
        # Independent ruler: brute-force minimum inter-molecular distance
        # under the minimum image must respect the tolerance.
        mol = np.repeat(np.arange(copies), n)
        d = pos[:, None, :] - pos[None, :, :]
        d -= np.round(d / edge) * edge
        dist = np.sqrt((d**2).sum(axis=-1))
        inter = mol[:, None] != mol[None, :]
        assert dist[inter].min() >= 2.0 - 1e-9

    def test_grow_pack_density_resolved_box(self):
        # ``with_density`` resolves the cubic periodic box from the frozen
        # stage-① formula, and ``Target.with_mass`` OVERRIDES the
        # element-derived mass (5 × C would be 60.055 amu, not 72).
        rho, mass_per_copy, copies = 0.05, 72.0, 2
        target = Target(_chain_frame(5), copies).with_mass(mass_per_copy)
        result = (
            _grow()
            .with_seed(3)
            .with_tolerance(2.0)
            .with_density(rho)
            .run([target], max_loops=50)
        )
        assert result.softened == 0
        assert result.fdist == 0.0

        expected_edge = (mass_per_copy * copies / (AVOGADRO * rho) * 1e24) ** (
            1.0 / 3.0
        )
        box = result.frame.box
        assert box is not None, "the density-resolved box must be stamped on the frame"
        np.testing.assert_allclose(
            np.asarray(box.lengths),
            [expected_edge] * 3,
            rtol=1e-9,
        )


class TestGencanPath:
    """The rigid-body entry through the same result type."""

    def test_gencan_softened_is_zero_and_deterministic(self):
        # A GENCAN pack reports softened == 0, and the same seed reproduces
        # the same positions bitwise — one entry per run, one verdict.
        def pack():
            return (
                GenCanPack()
                .with_seed(42)
                .with_tolerance(2.0)
                .with_progress(False)
                .with_periodic_box([0.0, 0.0, 0.0], [15.0, 15.0, 15.0])
                .run([Target(_chain_frame(2, bonds=False), 3)], max_loops=50)
            )

        a = pack()
        assert a.softened == 0
        b = pack()
        assert b.softened == 0
        np.testing.assert_array_equal(a.positions, b.positions)


class TestLatticeGrow:
    """Diamond-lattice growth entry (lattice-growth-phase spec)."""

    def test_lattice_bead_chain_constructive(self):
        n, copies, edge = 12, 8, 26.0
        result = (
            LatticeGrow(TorsionPrior.uniform())
            .with_seed(7)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [edge] * 3)
            .run([Target(_chain_frame(n), copies)], max_loops=60)
        )
        assert result.converged
        assert result.softened == 0
        assert result.fdist == 0.0
        assert result.positions.shape == (copies * n, 3)
        assert "LatticeGrow" in repr(LatticeGrow(TorsionPrior.uniform()))

    def test_lattice_requires_torsion_prior(self):
        with pytest.raises(TypeError):
            LatticeGrow()  # ty: ignore[missing-argument]

    def test_lattice_then_seeded_push_off(self):
        # Dense cell: the lattice generates where the continuum grinds; a
        # non-converged result chains into the explicit seeded push-off.
        n, copies, edge = 24, 20, 22.0
        target = lambda: Target(_chain_frame(n), copies)  # noqa: E731
        grown = (
            LatticeGrow(TorsionPrior.uniform())
            .with_seed(11)
            .with_tolerance(2.0)
            .with_periodic_box([0.0, 0.0, 0.0], [edge] * 3)
            .run([target()], max_loops=60)
        )
        assert grown.natoms == copies * n
        if not grown.converged:
            pushed = (
                GenCanPack()
                .seeded_from(grown)
                .with_seed(11)
                .with_tolerance(2.0)
                .run([target()], max_loops=120)
            )
            assert pushed.natoms == copies * n
            assert pushed.fdist <= grown.fdist
