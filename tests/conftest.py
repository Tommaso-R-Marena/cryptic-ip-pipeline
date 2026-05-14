import os
import sys
from pathlib import Path

# Ensure imports work without `pip install -e .`
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def tiny_pdb() -> Path:
    return FIXTURES / "tiny.pdb"


@pytest.fixture
def fpocket_info() -> Path:
    return FIXTURES / "fpocket_info.txt"


@pytest.fixture
def dx_file() -> Path:
    return FIXTURES / "potential.dx"


@pytest.fixture
def tmp_results(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "results").mkdir()
    return tmp_path
