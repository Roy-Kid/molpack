"""molpack version and molrs minor-line check.

The major.minor of the installed ``molcrafts-molrs`` must match the line this
wheel was built against. Patch may drift. A mismatch is an ``ImportError`` at
import time, not a later Frame-capsule segfault.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version

try:
    version = _pkg_version("molcrafts-molpack")
except PackageNotFoundError:
    version = "0.3.0"

# Keep in lockstep with ``python/pyproject.toml`` (``molcrafts-molrs>=X.Y,<X.Y+1``)
# and ``MOLRS_GIT_REF`` in ``.github/workflows/ci.yml``.
MOLRS_MINOR: tuple[int, int] = (0, 15)


def _minor_tuple(ver: str) -> tuple[int, int]:
    core = ver.split("+", 1)[0].split("-", 1)[0]
    parts = core.split(".")
    if len(parts) < 2:
        raise ValueError(f"expected at least major.minor, got {ver!r}")
    try:
        return int(parts[0]), int(parts[1])
    except ValueError as exc:
        raise ValueError(f"non-numeric version {ver!r}") from exc


def check_molrs_version() -> str:
    """Require the installed molrs to match ``MOLRS_MINOR``; return its version."""
    # Looked up at call time (not the module-level import) so tests can
    # monkeypatch ``importlib.metadata.version``.
    from importlib.metadata import version as pkg_version

    name = "molcrafts-molrs"
    try:
        installed = pkg_version(name)
    except PackageNotFoundError as exc:
        raise ImportError(
            f"molpack requires {name}, but its package metadata is missing"
        ) from exc

    try:
        got = _minor_tuple(installed)
    except ValueError as exc:
        raise ImportError(f"Cannot parse {name} version {installed!r} ({exc})") from exc

    if got == MOLRS_MINOR:
        return installed

    major, minor = MOLRS_MINOR
    raise ImportError(
        f"Minor-version mismatch: molpack expects {name} {major}.{minor}.* "
        f"but {name} {installed} is installed. Rebuild molpack against that "
        f"molrs (`maturin develop` in molpack/python) or install "
        f"`{name}>={major}.{minor}.0,<{major}.{minor + 1}`."
    )


check_molrs_version()
