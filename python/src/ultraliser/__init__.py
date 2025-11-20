"""Python interface for Ultraliser volume-to-mesh reconstruction."""

from __future__ import annotations

from importlib import metadata as _metadata

from ._core import volume_to_mesh

__all__ = ["volume_to_mesh"]

try:
    __version__ = _metadata.version("ultraliser-python")
except _metadata.PackageNotFoundError:  # pragma: no cover - during local development
    __version__ = "0.0.0"
