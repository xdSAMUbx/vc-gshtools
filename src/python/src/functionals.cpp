// functionals.cpp: the 45 geopotential functionals (functionals.hpp), height sources (height.hpp, height_model.hpp).
//
// Units of the returned values (SI, from functionals.hpp): V, W, T, dV_phi, dV_lambda in m^2/s^2 (angular derivatives
// per radian); dV_r, g_x, g_y, g_z, g, g_sa, dg_sa, Dg_sa, dg in m/s^2; second derivatives (V_*, T_*, W_rr) in 1/s^2
// (1 Eotvos = 1e-9 1/s^2); xi, eta, Theta in radians; zeta_ell, N, zeta in metres.
#include <cstdint>
#include <memory>
#include "common.hpp"
#include "functionals.hpp"
#include "height.hpp"
#include "height_model.hpp"
#include "model.hpp"
#include "synth_types.hpp"

namespace gshpy
{
  namespace
  {
// The 45 functionals, in enum order.
#define GSH_FUNCTIONAL_LIST(X)                                                                                      \
  X(V) X(dV_r) X(dV_phi) X(dV_lambda)                                                                               \
  X(W) X(g_x) X(g_y) X(g_z) X(g) X(g_sa) X(W_rr)                                                                    \
  X(V_rr) X(V_rphi) X(V_rlambda) X(V_phiphi) X(V_philambda) X(V_lambdalambda)                                       \
  X(V_xx) X(V_xy) X(V_xz) X(V_yy) X(V_yz) X(V_zz)                                                                   \
  X(T) X(dg_sa) X(Dg_sa) X(T_rr) X(dg)                                                                              \
  X(T_rphi) X(T_rlambda) X(T_phiphi) X(T_philambda) X(T_lambdalambda)                                               \
  X(T_xx) X(T_xy) X(T_xz) X(T_yy) X(T_yz) X(T_zz)                                                                   \
  X(xi) X(eta) X(Theta) X(zeta_ell) X(N) X(zeta)

    // compute_functionals result as seen from Python.
    struct PyFunctionalResult
    {
      py::dict values;       // {Functional: float64 array}
      std::int64_t points_without_height = 0;
      std::int64_t polar_points = 0;
      py::dict timing_s;     // {stage name: seconds}
    };
  }

  void bind_functionals(py::module_ &m)
  {
    using namespace gsh;

    // ---- Functional enum ---------------------------------------------------------------------------------------
    py::enum_<Functional> fe(m, "Functional", "The 45 geopotential functionals (Bucha and Janak 2019, GrafLab definitions).");
#define GSH_ENUM_VALUE(name) fe.value(#name, Functional::name);
    GSH_FUNCTIONAL_LIST(GSH_ENUM_VALUE)
#undef GSH_ENUM_VALUE
    m.def("functional_name", [](Functional f) { return std::string(functional_name(f)); }, py::arg("functional"),
          "Name of a functional as a string.");
    m.def("all_functionals", []() { return all_functionals(); }, "All 45 functionals, in enum order.");

    // ---- height sources ----------------------------------------------------------------------------------------
    py::class_<HeightSource, std::shared_ptr<HeightSource>>(m, "HeightSource",
                                                            "Abstract source of topographic heights H (metres).")
        .def("evaluate", [](const HeightSource &hs, const Points &pts, const DArr &lat_geodetic_deg, const SynthOptions &o) {
               auto la = copy_in<AlignedVec<double>>(lat_geodetic_deg, "lat_geodetic_deg");
               require_same(static_cast<ssize>(la.size()), static_cast<ssize>(pts.size()), "lat_geodetic_deg, pts");
               std::vector<double> H;
               {
                 py::gil_scoped_release g;
                 hs.evaluate(pts, la, H, o);
               }
               return to_numpy(std::move(H));
             },
             py::arg("pts"), py::arg("lat_geodetic_deg"), py::arg("opts") = SynthOptions(),
             "Heights [m] at the points; NaN where there is no data. lat_geodetic_deg: geodetic latitude in DEGREES.");

    py::class_<HeightGrid>(m, "HeightGrid", "Regular grid of heights (rows north to south, columns west to east).")
        .def(py::init<>())
        .def_readwrite("lat_south", &HeightGrid::lat_south)
        .def_readwrite("lat_north", &HeightGrid::lat_north)
        .def_readwrite("lon_west", &HeightGrid::lon_west)
        .def_readwrite("lon_east", &HeightGrid::lon_east)
        .def_readwrite("nlat", &HeightGrid::nlat)
        .def_readwrite("nlon", &HeightGrid::nlon)
        .def_property(
            "h", [](py::object self) { return view_of(self.cast<const HeightGrid &>().h, self); },
            [](HeightGrid &g, py::array_t<float, py::array::c_style | py::array::forcecast> a) {
              if (a.ndim() != 1)
                throw py::value_error("h must be one-dimensional");
              const float *p = a.data();
              g.h.assign(p, p + a.size());
            },
            "heights [m], float32, read-only view (assign a new array to replace)")
        .def("dlat", &HeightGrid::dlat, "grid step in latitude [deg]")
        .def("dlon", &HeightGrid::dlon, "grid step in longitude [deg]")
        .def("interpolate", &HeightGrid::interpolate, py::arg("lat_deg"), py::arg("lon_deg"),
             "bilinear height [m] at degrees; NaN outside the grid")
        .def("block_mean", &HeightGrid::block_mean, py::arg("block_deg"),
             "coarser grid of block averages over about block_deg degrees");

    m.def("load_gravsoft_grid", [](const std::string &path) {
           HeightGrid g;
           {
             py::gil_scoped_release gil;
             g = load_gravsoft_grid(path);
           }
           return g;
         },
          py::arg("path"), "Reads a GRAVSOFT grid (header lat1 lat2 lon1 lon2 dlat dlon). Raises RuntimeError on errors.");

    m.def("height_from_grid", [](HeightGrid grid, double block_deg) { return height_from_grid(std::move(grid), block_deg); },
          py::arg("grid"), py::arg("block_deg") = 0.0,
          "Height source that interpolates the grid. block_deg > 0 first block-averages it (degrees).");

    m.def("height_from_model", [](Model model, int nmax) { return height_from_model(std::move(model), nmax); },
          py::arg("model"), py::arg("nmax") = -1,
          "Height source from an HC/HS topography model (evaluated on the sphere, r = R).");

    // ---- FunctionalRequest and result --------------------------------------------------------------------------
    py::class_<FunctionalRequest>(m, "FunctionalRequest", "What compute_functionals computes and how. Units: see module docs.")
        .def(py::init<>())
        .def_readwrite("outputs", &FunctionalRequest::outputs, "list of Functional to evaluate")
        .def_readwrite("nmin", &FunctionalRequest::nmin, "lowest degree included")
        .def_readwrite("nmax", &FunctionalRequest::nmax, "highest degree; < 0 = model degree")
        .def_readwrite("ellipsoid", &FunctionalRequest::ellipsoid, "reference ellipsoid (default GRS80)")
        .def_readwrite("height", &FunctionalRequest::height, "HeightSource; required for N and zeta, else None")
        .def_readwrite("density", &FunctionalRequest::density, "topographic density [kg/m^3]")
        .def_readwrite("G", &FunctionalRequest::G, "gravitational constant [m^3 kg^-1 s^-2]")
        .def_readwrite("pole_cos_threshold", &FunctionalRequest::pole_cos_threshold,
                       "below this cos(phi) the oracle is used")
        .def_readwrite("synth", &FunctionalRequest::synth, "SynthOptions (threads, ISA, switches)");

    py::class_<PyFunctionalResult>(m, "FunctionalResult", "Result of compute_functionals.")
        .def_readonly("values", &PyFunctionalResult::values, "dict {Functional: float64 array}, one entry per point")
        .def_readonly("points_without_height", &PyFunctionalResult::points_without_height,
                      "points where the height source returned NaN (N and zeta are NaN there)")
        .def_readonly("polar_points", &PyFunctionalResult::polar_points, "points evaluated with the oracle")
        .def_readonly("timing_s", &PyFunctionalResult::timing_s, "dict {stage: wall seconds}");

    m.def("compute_functionals", [](const Model &model, const Points &pts, const FunctionalRequest &req) {
           FunctionalResult res;
           {
             py::gil_scoped_release g;
             res = gsh::compute_functionals(model, pts, req);
           }
           PyFunctionalResult out;
           for (auto &kv : res.values)
             out.values[py::cast(kv.first)] = to_numpy(std::move(kv.second));
           out.points_without_height = res.points_without_height;
           out.polar_points = res.polar_points;
           for (const auto &kv : res.timing_s)
             out.timing_s[py::str(kv.first)] = kv.second;
           return out;
         },
          py::arg("model"), py::arg("pts"), py::arg("req"),
          "Computes the requested functionals at every point, with the fewest synthesis passes.\n"
          "Raises ValueError: nmax above model.nmax; N or zeta without height; a functional undefined at a pole.");
  }
}
