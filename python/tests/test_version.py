"""molpack/molrs minor-version compatibility tests."""

from __future__ import annotations

import importlib
import importlib.metadata

import pytest

version_module = importlib.import_module("molpack.version")


def test_matching_minor_is_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    major, minor = version_module.MOLRS_MINOR
    molrs_ver = f"{major}.{minor}.0"

    def fake_version(name: str) -> str:
        if name == "molcrafts-molrs":
            return molrs_ver
        raise importlib.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(importlib.metadata, "version", fake_version)
    assert version_module.check_molrs_version() == molrs_ver


def test_same_minor_different_patch_is_accepted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    major, minor = version_module.MOLRS_MINOR
    for patch in ("0", "1", "99"):
        molrs_ver = f"{major}.{minor}.{patch}"
        monkeypatch.setattr(importlib.metadata, "version", lambda _n, v=molrs_ver: v)
        assert version_module.check_molrs_version() == molrs_ver


def test_different_minor_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(importlib.metadata, "version", lambda _name: "0.12.0")
    with pytest.raises(ImportError, match="Minor-version mismatch"):
        version_module.check_molrs_version()


def test_missing_molrs_metadata_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    def missing(_name: str) -> str:
        raise importlib.metadata.PackageNotFoundError("molcrafts-molrs")

    monkeypatch.setattr(importlib.metadata, "version", missing)
    with pytest.raises(ImportError, match="package metadata is missing"):
        version_module.check_molrs_version()
