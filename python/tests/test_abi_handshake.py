"""FFI ABI gates: versioned capsule names + the import-time handshake.

Minor-line = ABI version: molpack exchanges ``molrs_ffi`` handle capsules with
the installed ``molcrafts-molrs`` wheel, so both must embed the same molrs
``major.minor``. The import-time handshake (``interop::check_abi``) already
passed — this module is importable — so these tests pin the *other* gate: a
capsule from a different minor line (spelled with the pre-0.14 unversioned
name here) must be rejected at resolve time, cleanly.
"""

from __future__ import annotations

import ctypes

import molrs
import pytest

import molpack


class _LegacyFrame:
    """Quacks like a molrs Frame but exports a pre-0.14 unversioned capsule."""

    def _ffi_frameref_capsule(self):  # noqa: ANN202 — mirrors the duck-typed contract
        new_capsule = ctypes.pythonapi.PyCapsule_New
        new_capsule.restype = ctypes.py_object
        new_capsule.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
        # Bogus payload on purpose: the name check must reject the capsule
        # before any dereference happens.
        return new_capsule(ctypes.c_void_p(0xDEAD), b"molrs.FrameRef", None)


class TestVersionedCapsuleGate:
    def test_cross_minor_capsule_is_rejected_at_resolve(self) -> None:
        with pytest.raises(ValueError, match="minor line"):
            molpack.Target(_LegacyFrame(), 1)

    def test_same_line_frame_resolves(self) -> None:
        import numpy as np

        frame = molrs.Frame()
        block = molrs.Block()
        block.insert("x", np.array([0.0, 1.0]))
        block.insert("y", np.array([0.0, 0.0]))
        block.insert("z", np.array([0.0, 0.0]))
        block.insert("id", np.array([1, 2], dtype=np.uint32))
        frame["atoms"] = block
        assert molpack.Target(frame, 1) is not None
