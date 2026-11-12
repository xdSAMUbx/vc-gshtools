"""gshtools: Python bindings of the GSH_Tools C++ library.

The C++ library is the reference implementation. This package only converts types and calls it; no numerical
logic is implemented here. Units are documented on each function (angles in radians unless the name says _deg).
"""

import os as _os
import sys as _sys

__version__ = "2.0.0"

# Windows: the extension needs the DLLs that were installed next to it (the OpenMP runtime, see
# python/CMakeLists.txt). Python 3.8+ does not search the package directory for them, so it is added explicitly.
_dll_dir_handle = None
if _sys.platform == "win32" and hasattr(_os, "add_dll_directory"):
    _dll_dir_handle = _os.add_dll_directory(_os.path.dirname(_os.path.abspath(__file__)))

from ._core import *  # noqa: E402,F401,F403  (re-exports every public name of the extension)
from . import _core as _core  # noqa: E402

# Submodules are created by pybind11 as attributes; register them so that `import gshtools.tides` also works.
for _name in ("coords", "tides", "legendre"):
    _sys.modules[f"{__name__}.{_name}"] = getattr(_core, _name)
del _name
