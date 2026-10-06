"""The Python mirror of the multi-stage ``Pipeline`` (stage-pipeline-07-bindings).

``Pipeline`` is the composition surface: the entries (``GenCanPack`` /
``CbmcGrow`` / ``LatticeGrow``) stay single-stage presets, and chaining them
is one object with the *same* shared knobs. What this file owns, and nothing
else can:

1. **Two spellings, one answer.** A single-stage pipeline is the preset run,
   bitwise (`ac-001`); a two-stage pipeline runs and reports one verdict.
2. **Nothing is dropped silently.** A stage's own ``with_handler`` callback is
   *adopted* by the pipeline and still fires (`ac-004`); a preset carrying a
   non-default *shared* knob into a pipeline is refused by name, knob included
   (`ac-004`); an empty pipeline is a named ``ValueError``, not a no-op; an
   object that is not a registered entry is a ``TypeError`` that *lists* the
   entries (`ac-002`'s user-visible face).
3. **Handlers can see which stage they are in.** ``StepInfo.stage`` carries the
   ``index`` / ``total`` / ``name`` triple, in a pipeline and in a bare preset
   run alike (`ac-003`).
4. **The numbers do not drift.** Two hard-coded goldens (`ac-007`).

Fixtures are copied from ``test_packer.py`` / ``test_grow.py`` rather than
shared, so a change to either file cannot silently move this file's answers.
Deterministic by construction: fixed seeds, no wall clock, no filesystem, no
network, no third-party oracle.

Gate: ``uv run --directory python --group dev tox -e py``.
"""

from __future__ import annotations

import molrs
import numpy as np
import pytest

from molpack import (
    CbmcGrow,
    GenCanPack,
    LatticeGrow,
    Pipeline,
    StepContext,
    StepInfo,
    Target,
    TorsionPrior,
)

#: The entries a `Pipeline` accepts as a stage. The binding builds its
#: `TypeError` text from ONE registry; this tuple is the same list spelled
#: from the classes themselves, so it cannot drift into a stale literal.
ENTRY_NAMES = (GenCanPack.__name__, CbmcGrow.__name__, LatticeGrow.__name__)


# ── fixtures ──────────────────────────────────────────────────────────────


def _water_frame() -> molrs.Frame:
    """Rigid water template — the exact literals of ``tests/pipeline.rs``'s
    ``water()`` fixture, so the Rust and Python goldens describe one run."""
    return molrs.Frame(
        {
            "atoms": {
                "x": np.array([0.0, 0.96, -0.24], dtype=np.float64),
                "y": np.array([0.0, 0.0, 0.93], dtype=np.float64),
                "z": np.zeros(3, dtype=np.float64),
                "element": ["O", "H", "H"],
            }
        }
    )


def _water_targets(count: int = 60) -> list[Target]:
    """``count`` waters held by a plain (non-periodic) box restraint — no box
    and no cell declaration, so the packing volume is the fall-back one."""
    return [
        Target(_water_frame(), count)
        .with_name("water")
        .with_restraint(molrs.Cuboid([0.0, 0.0, 0.0], [14.0, 14.0, 14.0]))
    ]


#: A planar zigzag bead chain written as *literals*, never as
#: ``bond * cos(alpha)``: a golden must not depend on which libm rounded the
#: fixture. Bond lengths are ~1.5 Å by construction (1.2² + 0.9² = 1.5²).
_CHAIN_X = (0.0, 1.2, 2.4, 3.6, 4.8)
_CHAIN_Z = (0.0, 0.9, 0.0, 0.9, 0.0)


def _chain_frame() -> molrs.Frame:
    """Bonded 5-bead chain — the smallest growable fixture ``test_grow.py``
    uses, with the coordinates pinned to literals."""
    n = len(_CHAIN_X)
    return molrs.Frame(
        {
            "atoms": {
                "x": np.array(_CHAIN_X, dtype=np.float64),
                "y": np.zeros(n, dtype=np.float64),
                "z": np.array(_CHAIN_Z, dtype=np.float64),
                "element": ["C"] * n,
            },
            "bonds": {
                "atomi": np.arange(0, n - 1, dtype=np.uint64),
                "atomj": np.arange(1, n, dtype=np.uint64),
            },
        }
    )


def _chain_targets(copies: int = 2) -> list[Target]:
    return [Target(_chain_frame(), copies)]


class _Counter:
    """Counting handler (the shape ``test_packer.py`` / ``test_handler.py``
    use): three optional hooks, no state beyond the tallies."""

    def __init__(self) -> None:
        self.starts = 0
        self.steps = 0
        self.finishes = 0

    def on_start(self, ntotat: int, ntotmol: int) -> None:
        self.starts += 1

    def on_step(self, info: StepInfo, ctx: StepContext) -> None:
        self.steps += 1

    def on_finish(self) -> None:
        self.finishes += 1


class _StageRecorder:
    """Records the ``StepInfo.stage`` triple of every step it sees."""

    def __init__(self) -> None:
        self.seen: list[tuple[int, int, str]] = []

    def on_step(self, info: StepInfo, ctx: StepContext) -> None:
        self.seen.append((info.stage.index, info.stage.total, info.stage.name))


# ── happy path ────────────────────────────────────────────────────────────


def test_pipeline_two_stages_runs() -> None:
    """Growth then rigid push-off is ONE object with one verdict."""
    result = (
        Pipeline([CbmcGrow(TorsionPrior.uniform()), GenCanPack()])
        .with_seed(7)
        .with_tolerance(1.0)
        .with_periodic_box([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])
        .run(_chain_targets(), max_loops=50)
    )

    assert result.natoms == 10  # 2 copies × 5 beads
    assert result.positions.shape == (10, 3)
    assert np.isfinite(result.positions).all()
    assert isinstance(result.converged, bool)


def test_pipeline_adopts_stage_handlers() -> None:
    """A handler that arrived on a stage observes the whole run.

    Handing a preset to a pipeline must not cost the caller their callback:
    the stage's handlers are adopted, bracketed once by the run.
    """
    counter = _Counter()
    Pipeline([GenCanPack().with_handler(counter)]).with_seed(42).with_tolerance(
        2.0
    ).run(_water_targets(), max_loops=20)

    assert counter.steps > 0, (
        "the handler carried in by GenCanPack().with_handler saw no on_step "
        "events — a pipeline adopts a stage's handlers, it does not drop them"
    )
    assert counter.starts == 1
    assert counter.finishes == 1


def test_step_info_exposes_stage_triple() -> None:
    """``StepInfo.stage`` names the emitting stage, in both spellings."""
    recorder = _StageRecorder()
    Pipeline([CbmcGrow(TorsionPrior.uniform()), GenCanPack()]).with_handler(
        recorder
    ).with_seed(7).with_tolerance(1.0).with_periodic_box(
        [0.0, 0.0, 0.0], [20.0, 20.0, 20.0]
    ).run(_chain_targets(), max_loops=50)

    assert recorder.seen, "expected at least one on_step call"
    assert {total for _, total, _ in recorder.seen} == {2}, (
        "every step of a two-stage run reports total == 2"
    )
    indices = [index for index, _, _ in recorder.seen]
    assert indices == sorted(indices), f"stage index went backwards: {indices}"
    assert set(indices) <= {0, 1}
    assert {name for _, _, name in recorder.seen} <= {"growth", "gencan"}

    # A bare preset run is a one-stage run and says so.
    solo = _StageRecorder()
    GenCanPack().with_handler(solo).with_seed(42).with_tolerance(2.0).run(
        _water_targets(), max_loops=20
    )
    assert solo.seen, "expected at least one on_step call"
    assert set(solo.seen) == {(0, 1, "gencan")}


# ── named refusals ────────────────────────────────────────────────────────


def test_pipeline_preset_settings_inside_pipeline_is_value_error() -> None:
    """Shared knobs are the run's, and a refusal names the stage and knob."""
    with pytest.raises(ValueError) as excinfo:
        Pipeline([GenCanPack().with_seed(7)]).run(_water_targets(4), max_loops=2)

    message = str(excinfo.value)
    assert "gencan" in message, f"the refusal must name the stage, got: {message}"
    assert "seed" in message, f"the refusal must name the knob, got: {message}"

    # The supported spelling of the same intent: the knob on the Pipeline.
    result = (
        Pipeline([GenCanPack()])
        .with_seed(7)
        .with_tolerance(2.0)
        .run(_water_targets(4), max_loops=2)
    )
    assert result.natoms == 12


def test_pipeline_empty_is_value_error() -> None:
    """An empty pipeline is a named error, never a silent no-op run.

    This is the composition error Python *can* reach. The sibling
    ``PackError::StageOrder`` ("stage X requires placements but nothing before
    it placed the molecules") is unreachable from Python by construction: the
    only stages Python can build are the three registered entries, and every
    one of them places molecules itself (``Requires::nothing``). Building a
    stage that *requires* prior placements is a Rust-level extension point
    (`Stage` stays Rust-only, spec §Design), so `StageOrder` is owned by
    ``tests/pipeline.rs``.
    """
    for empty in (Pipeline([]), Pipeline()):
        with pytest.raises(ValueError, match="no stages"):
            empty.run(_water_targets(2), max_loops=1)


def test_pipeline_rejects_unknown_stage_with_registry_message() -> None:
    """A non-entry object is a ``TypeError`` that lists the entries.

    The list comes from the binding's one stage registry, so a fourth entry
    shows up here without anyone editing an error string.
    """
    # Deliberate misuse: the stub's `StageEntry` union is closed on purpose, so
    # the checker is right and the runtime TypeError is the contract under test.
    with pytest.raises(TypeError) as from_ctor:
        Pipeline([object()])  # ty: ignore[invalid-argument-type]

    with pytest.raises(TypeError) as from_builder:
        Pipeline().with_stage(42)  # ty: ignore[invalid-argument-type]

    for excinfo in (from_ctor, from_builder):
        message = str(excinfo.value)
        for name in ENTRY_NAMES:
            assert name in message, (
                f"the refusal must list the supported entries; {name!r} missing "
                f"from: {message}"
            )
