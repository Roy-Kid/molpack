"""Tests for Python-defined packing handlers (``Molpack.with_handler``)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest

import molpack


def _two_water_frame() -> dict:
    """Two trivially distinct atoms so packing has something to do."""
    return {
        "atoms": {
            "x": np.array([0.0, 1.5]),
            "y": np.array([0.0, 0.0]),
            "z": np.array([0.0, 0.0]),
            "element": ["O", "H"],
        }
    }


def _packer() -> molpack.Molpack:
    return molpack.Molpack().with_progress(False).with_inner_iterations(5)


@dataclass
class CallLog:
    started: bool = False
    finished: bool = False
    ntotat: int = 0
    ntotmol: int = 0
    steps: list = field(default_factory=list)

    def on_start(self, ntotat: int, ntotmol: int) -> None:
        self.started = True
        self.ntotat = ntotat
        self.ntotmol = ntotmol

    def on_step(self, info: molpack.StepInfo) -> None:
        self.steps.append((info.phase, info.loop_idx, info.fdist, info.frest))

    def on_finish(self) -> None:
        self.finished = True


class TestHandlerCallbacks:
    def test_lifecycle_methods_are_called(self):
        log = CallLog()
        # 30 copies in a 6³ box forces the outer loop to actually
        # iterate — a too-easy setup converges during init and
        # `on_step` never fires.
        target = molpack.Target(_two_water_frame(), count=30).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [6.0, 6.0, 6.0])
        )
        _packer().with_handler(log).with_seed(1).pack([target], max_loops=5)

        assert log.started is True
        assert log.finished is True
        # 30 copies × 2 atoms per copy → 60 atoms; 30 molecules.
        assert log.ntotat == 60
        assert log.ntotmol == 30
        assert len(log.steps) >= 1

    def test_step_info_fields_are_readable(self):
        captured: list[tuple] = []

        class Grabber:
            def on_step(self, info: molpack.StepInfo) -> None:
                captured.append(
                    (
                        info.loop_idx,
                        info.max_loops,
                        info.phase,
                        info.total_phases,
                        info.molecule_type,
                        info.fdist,
                        info.frest,
                        info.improvement_pct,
                        info.radscale,
                        info.precision,
                        info.relaxer_acceptance,
                    )
                )

        target = molpack.Target(_two_water_frame(), count=2).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [5.0, 5.0, 5.0])
        )
        _packer().with_handler(Grabber()).with_seed(1).pack([target], max_loops=2)

        assert captured, "expected at least one on_step call"
        first = captured[0]
        assert isinstance(first[0], int)  # loop_idx
        assert first[1] == 2  # max_loops
        assert first[2] >= 0  # phase
        assert first[3] >= 1  # total_phases
        # molecule_type is None for the final all-types phase, int otherwise
        assert first[4] is None or isinstance(first[4], int)
        assert all(isinstance(v, float) for v in first[5:10])
        # relaxer_acceptance is a (possibly empty) list of (int, float) tuples
        assert isinstance(first[10], list)

    def test_methods_are_optional(self):
        """Handler with *no* methods must not crash."""

        class Empty:
            pass

        target = molpack.Target(_two_water_frame(), count=2).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [5.0, 5.0, 5.0])
        )
        result = (
            _packer().with_handler(Empty()).with_seed(1).pack([target], max_loops=2)
        )
        assert result.natoms == 4


class TestHandlerEarlyStop:
    def test_returning_true_halts_pack(self):
        steps_seen: list[int] = []

        class StopAfterOne:
            def on_step(self, info: molpack.StepInfo) -> bool:
                steps_seen.append(info.loop_idx)
                return True  # request immediate stop

        target = molpack.Target(_two_water_frame(), count=4).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
        )
        _packer().with_handler(StopAfterOne()).with_seed(1).pack([target], max_loops=50)

        # The per-phase compaction loop itself runs through its handler
        # pass before checking should_stop; we just assert that we did
        # NOT run the full 50 outer loops.
        assert len(steps_seen) < 50, (
            f"early stop was ignored; saw {len(steps_seen)} steps"
        )


class TestHandlerErrorPropagation:
    def test_exception_in_on_step_is_reraised(self):
        class Explodes:
            def on_step(self, info: molpack.StepInfo) -> None:
                raise ValueError("boom from handler")

        target = molpack.Target(_two_water_frame(), count=2).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [5.0, 5.0, 5.0])
        )
        with pytest.raises(ValueError, match="boom from handler"):
            _packer().with_handler(Explodes()).with_seed(1).pack([target], max_loops=5)

    def test_exception_in_on_start_is_reraised(self):
        class ExplodesEarly:
            def on_start(self, ntotat: int, ntotmol: int) -> None:
                raise RuntimeError("boom from on_start")

        target = molpack.Target(_two_water_frame(), count=2).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [5.0, 5.0, 5.0])
        )
        with pytest.raises(RuntimeError, match="boom from on_start"):
            _packer().with_handler(ExplodesEarly()).with_seed(1).pack(
                [target], max_loops=5
            )


class TestMultipleHandlers:
    def test_each_handler_receives_callbacks(self):
        log1 = CallLog()
        log2 = CallLog()
        target = molpack.Target(_two_water_frame(), count=2).with_restraint(
            molpack.InsideBoxRestraint([0.0, 0.0, 0.0], [5.0, 5.0, 5.0])
        )
        _packer().with_handler(log1).with_handler(log2).with_seed(1).pack(
            [target], max_loops=2
        )

        assert log1.started and log2.started
        assert log1.finished and log2.finished
        assert len(log1.steps) == len(log2.steps) >= 1
