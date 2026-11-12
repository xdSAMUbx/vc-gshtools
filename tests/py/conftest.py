"""Shared fixtures for the gshtools tests.

The tests need the compiled extension. If `import gshtools` fails (in-tree build, not installed), the in-tree
build folder build/py/python is added to sys.path. Tests that need the EGM2008 model are skipped when the file is
missing.
"""

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
EGM2008 = ROOT / "data" / "cache" / "modelos" / "EGM2008" / "EGM2008.gfc"
INTREE = ROOT / "build" / "py" / "python"

try:
    import gshtools  # noqa: F401
except ImportError:
    sys.path.insert(0, str(INTREE))


@pytest.fixture(scope="session")
def gsh():
    import gshtools

    return gshtools


@pytest.fixture(scope="session")
def egm2008_nmax60(gsh):
    """EGM2008 truncated at degree 60 (fast enough for the test suite)."""
    if not EGM2008.exists():
        pytest.skip(f"model file missing: {EGM2008}")
    return gsh.load_gfc(str(EGM2008), 60)


@pytest.fixture(scope="session")
def egm2008_sigma_nmax60(gsh):
    if not EGM2008.exists():
        pytest.skip(f"model file missing: {EGM2008}")
    return gsh.load_gfc_sigma(str(EGM2008), 60)
