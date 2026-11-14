// legendre.cpp: submodule legendre (legendre.hpp). Fully normalised associated Legendre functions.
//
// Arguments of sfc/msfc/era/compute: t = sin(phi) and u = cos(phi), phi = GEOCENTRIC latitude (radians).
// The result is a packed triangle with tri_size(nmax) entries; element (n, m) is at idx_tri(n, m, nmax).
#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <vector>
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

namespace gshpy
{
  namespace
  {
    // Python-side result of legendre_points(). The members are NumPy arrays, owned by this object.
    struct LegendrePoints
    {
      py::object P, dP, n, m;
    };

    // Moves a vector to NumPy as a C-contiguous (rows, cols) array without copying its elements.
    template <class V>
    py::array to_numpy_2d(V v, ssize rows, ssize cols)
    {
      using T = typename V::value_type;
      V *holder = new V(std::move(v));
      py::capsule owner(holder, [](void *p) { delete static_cast<V *>(p); });
      return py::array_t<T>(std::vector<ssize>{rows, cols},
                            std::vector<ssize>{cols * static_cast<ssize>(sizeof(T)), static_cast<ssize>(sizeof(T))},
                            holder->data(), owner);
    }
  }

  void bind_legendre_points(py::module_ &mod)
  {
    namespace L = gsh::legendre;

    py::class_<LegendrePoints>(mod, "LegendrePoints",
                               "Result of legendre_points(): P and dP are (npoints, count) float64 arrays (dP is None unless "
                               "derivative=True); n and m are int32 arrays of length count, the degree and order of each column.")
        .def_readonly("P", &LegendrePoints::P)
        .def_readonly("dP", &LegendrePoints::dP)
        .def_readonly("n", &LegendrePoints::n)
        .def_readonly("m", &LegendrePoints::m);

    mod.def(
        "legendre_points",
        [](const DArr &lat, int nmax, std::optional<std::vector<int>> orders, int nmin, bool derivative, int threads,
           py::object coefs) {
          const auto latv = copy_in<std::vector<double>>(lat, "lat");
          L::Selection sel;
          sel.orders = orders.value_or(std::vector<int>{});
          sel.nmin = nmin;
          sel.nmax = nmax;
          sel.derivative = derivative;
          // A user-supplied Coefs is used as is; py::object keeps it alive while the GIL is released.
          const bool given = !coefs.is_none();
          const L::Coefs *user = given ? &coefs.cast<const L::Coefs &>() : nullptr;
          L::Extracted ex;
          {
            py::gil_scoped_release g;
            if (given)
              ex = L::extract(std::span<const double>(latv.data(), latv.size()), *user, sel, threads);
            else
            {
              const L::Coefs c = L::make_coefs(nmax);
              ex = L::extract(std::span<const double>(latv.data(), latv.size()), c, sel, threads);
            }
          }
          const auto rows = static_cast<ssize>(latv.size());
          const auto cols = static_cast<ssize>(ex.count);
          LegendrePoints r;
          r.P = to_numpy_2d(std::move(ex.P), rows, cols);
          r.dP = derivative ? to_numpy_2d(std::move(ex.dP), rows, cols) : py::object(py::none());
          r.n = to_numpy(std::move(ex.n));
          r.m = to_numpy(std::move(ex.m));
          return r;
        },
        py::arg("lat"), py::arg("nmax"), py::arg("orders") = py::none(), py::arg("nmin") = 0,
        py::arg("derivative") = false, py::arg("threads") = 0, py::arg("coefs") = py::none(),
        "Fully normalised Legendre functions for many points, keeping only the selected orders and degrees.\n\n"
        "lat: geocentric latitudes in radians (array or sequence). nmax: maximum degree (<= coefs.nmax if coefs is given).\n"
        "orders: list of orders m (None = all 0..nmax); nmin: smallest degree (for order m the degrees stored are "
        "max(m, nmin)..nmax).\n"
        "derivative: also return dP/dphi. threads: 0 = OpenMP default. coefs: optional result of make_coefs(nmax) to reuse.\n"
        "Each point's full triangle is computed per thread and only the selected entries are copied out. Columns are in "
        "ascending m, then ascending n. Memory: 8 * npoints * count bytes (count = number of columns, times 2 with "
        "derivative). Returns a LegendrePoints object.");
  }
}
