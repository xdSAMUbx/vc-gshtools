// module.cpp: entry point of the extension module gshtools._core. Registration order matters only for default
// arguments and for types used as defaults; the order below respects that.
#include "common.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
  m.doc() = "Bindings of the GSH_Tools C++20 library. The C++ library is the reference implementation; "
            "this module only converts types.";

  gshpy::bind_model_points(m);
  gshpy::bind_ellipsoid(m);

  py::module_ coords = m.def_submodule("coords", "Coordinate conversions (degrees, metres).");
  gshpy::bind_coords(coords);

  py::module_ tides = m.def_submodule("tides", "Permanent-tide systems: C20, geoid height, height and gravity shifts.");
  gshpy::bind_tides(tides);

  py::module_ legendre = m.def_submodule("legendre", "Fully normalised associated Legendre functions (radians).");
  gshpy::bind_legendre(legendre);
  gshpy::bind_legendre_points(m);

  gshpy::bind_functionals(m);
  gshpy::bind_synthesis(m);
  gshpy::bind_errors(m);
}
