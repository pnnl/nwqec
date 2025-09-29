"""Python package for the NWQEC quantum transpiler bindings."""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("nwqec")
except _metadata.PackageNotFoundError:  # pragma: no cover - happens in local dev
    __version__ = "0.1.0-dev"

from ._core import *  # type: ignore[F401,F403]

__all__ = [name for name in globals() if not name.startswith("_")]
