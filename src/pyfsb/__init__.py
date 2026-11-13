"""pyfsb: a pymaster extension for Filtered-Squared Bispectra (FSBs)."""

from importlib.metadata import PackageNotFoundError, version

from .fsb import FSB, get_filters

try:
    __version__ = version("pyfsb")
except PackageNotFoundError:  # running from a source tree without installing
    __version__ = "unknown"

__all__ = ["FSB", "get_filters"]
