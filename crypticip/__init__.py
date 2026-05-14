"""crypticip — production pipeline for cryptic IP binding-site detection.

The public surface is intentionally small: most users interact through the
``crypticip`` CLI (see ``crypticip.cli``). Module-level imports here are kept
light so ``import crypticip`` works without optional binaries installed.
"""

__version__ = "0.5.0"

from .paths import ProjectPaths  # noqa: F401
from .config import load_config, Config  # noqa: F401
