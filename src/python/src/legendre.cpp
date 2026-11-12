// legendre.cpp: submodule legendre (legendre.hpp). Fully normalised associated Legendre functions.
//
// Arguments of sfc/msfc/era/compute: t = sin(phi) and u = cos(phi), phi = GEOCENTRIC latitude (radians).
// The result is a packed triangle with tri_size(nmax) entries; element (n, m) is at idx_tri(n, m, nmax).
#include <cstddef>
#include <span>
#include "common.hpp"
#include "legendre.hpp"

namespace gshpy
{
  void bind_legendre(py::module_ &m)
  {
    using namespace gsh;
    namespace L = gsh::legendre;

    py::class_<L::Coefs>(m, "Coefs", "Recurrence coefficients for one maximum degree (from make_coefs).")
        .def_readonly("nmax", &L::Coefs::nmax, "maximum degree");

    m.def("make_coefs", [](int nmax) {
           L::Coefs c;
           {
             py::gil_scoped_release g;
             c = L::make_coefs(nmax);
           }
           return c;
         },
          py::arg("nmax"), "Builds the recurrence coefficients for degree nmax. Cost O(nmax^2).");

    // The four evaluators share one shape: (t, u, coefs) -> packed triangle as a new array.
    using Evaluator = void (*)(double, double, const L::Coefs &, std::span<double>);
    auto def_eval = [&m](const char *name, Evaluator fn, const char *doc) {
      m.def(name, [fn](double t, double u, const L::Coefs &c) {
             AlignedVec<double> P(static_cast<std::size_t>(tri_size(c.nmax)));
             {
               py::gil_scoped_release g;
               fn(t, u, c, std::span<double>(P.data(), P.size()));
             }
             return to_numpy(std::move(P));
           },
           py::arg("t"), py::arg("u"), py::arg("coefs"), doc);
    };
    def_eval("sfc", &L::sfc, "Standard forward column (valid to degree 1800). t = sin(phi), u = cos(phi).");
    def_eval("msfc", &L::msfc, "Modified standard forward column (any degree; mantissa with R^-k scaling).");
    def_eval("era", &L::era, "Extended-range arithmetic reference (slow; precision reference).");
    def_eval("compute", &L::compute, "Dispatcher: sfc for nmax <= 1800, msfc otherwise.");

    m.def("dphi", [](int nmax, const DArr &P) {
           if (P.size() < tri_size(nmax))
             throw py::value_error("P must hold at least tri_size(nmax) elements");
           AlignedVec<double> dP(static_cast<std::size_t>(tri_size(nmax)));
           {
             py::gil_scoped_release g;
             L::dphi(nmax, std::span<const double>(P.data(), static_cast<std::size_t>(P.size())),
                     std::span<double>(dP.data(), dP.size()));
           }
           return to_numpy(std::move(dP));
         },
          py::arg("nmax"), py::arg("P"),
          "Latitude derivative dP/dphi of a packed triangle P (radians), non-singular at the poles. Returns a new array.");
  }
}
