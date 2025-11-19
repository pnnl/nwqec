"""Python package for the NWQEC quantum transpiler bindings."""

from importlib import metadata as _metadata

try:
    __version__ = _metadata.version("nwqec")
except _metadata.PackageNotFoundError:  # pragma: no cover - happens in local dev
    __version__ = "0.1.0-dev"

from ._core import *  # type: ignore[F401,F403]

# Platform-specific gridsynth availability
import platform as _platform
_is_windows = _platform.system() == "Windows"

# Check if C++ gridsynth is available
try:
    from ._core import WITH_GRIDSYNTH_CPP
    _has_cpp_gridsynth = WITH_GRIDSYNTH_CPP
except ImportError:
    _has_cpp_gridsynth = False

# Provide fallback message for Windows users
if _is_windows and not _has_cpp_gridsynth:
    import warnings
    warnings.warn(
        "C++ gridsynth backend not available on Windows. "
        "For gridsynth functionality, please install pygridsynth: pip install pygridsynth",
        UserWarning,
        stacklevel=2
    )

__all__ = [name for name in globals() if not name.startswith("_")]
