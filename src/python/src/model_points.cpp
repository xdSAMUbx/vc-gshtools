// model_points.cpp: Model and its reader, Points, SynthOptions, Isa, and read_table (table_io.hpp).
#include "common.hpp"
#include "ellipsoid.hpp"
#include "model.hpp"
#include "synth_types.hpp"
#include "table_io.hpp"

namespace gshpy
{
  namespace
  {
    // Table as seen from Python: columns is a list of float64 arrays, skipped_lines a list of 1-based line numbers.
    struct PyTable
    {
      py::list columns;
      std::vector<std::int64_t> skipped_lines;
    };

    void check_index(const gsh::Model &mo, int n, int m)
    {
      if (n < 0 || m < 0 || m > n || n > mo.nmax)
        throw py::index_error("require 0 <= m <= n <= nmax (model nmax = " + std::to_string(mo.nmax) + ")");
    }
  }

  void bind_model_points(py::module_ &m)
  {
    using namespace gsh;

    // ---- Model -------------------------------------------------------------------------------------------------
    py::class_<Model>(m, "Model", "Spherical-harmonic coefficient set (fully normalised). GM [m^3/s^2], R [m].")
        .def(py::init<>())
        .def_readonly("GM", &Model::GM, "gravitational constant [m^3/s^2]")
        .def_readonly("R", &Model::R, "reference radius [m]")
        .def_readonly("nmax", &Model::nmax, "maximum degree stored")
        .def_readonly("file_nmax", &Model::file_nmax, "max_degree declared in the file header")
        .def_readonly("name", &Model::name, "header modelname (may be empty)")
        .def_readonly("norm", &Model::norm, "header norm (may be empty)")
        .def_readonly("tide_system", &Model::tide_system, "header tide_system (may be empty)")
        .def_property_readonly("C", [](py::object self) { return view_of(self.cast<const Model &>().C, self); },
                               "packed C-bar (read-only view), index idx_tri(n, m, nmax)")
        .def_property_readonly("S", [](py::object self) { return view_of(self.cast<const Model &>().S, self); },
                               "packed S-bar (read-only view), index idx_tri(n, m, nmax)")
        .def("c", [](const Model &mo, int n, int k) { check_index(mo, n, k); return mo.c(n, k); }, py::arg("n"), py::arg("m"),
             "C-bar(n, m) [dimensionless]")
        .def("s", [](const Model &mo, int n, int k) { check_index(mo, n, k); return mo.s(n, k); }, py::arg("n"), py::arg("m"),
             "S-bar(n, m) [dimensionless]")
        .def("truncated", &Model::truncated, py::arg("new_nmax"), "copy with maximum degree reduced to new_nmax")
        .def_static("zeros", &Model::zeros, py::arg("GM"), py::arg("R"), py::arg("nmax"), "model with all coefficients zero");

    m.def("load_gfc", [](const std::string &path, int nmax) {
           py::gil_scoped_release g;
           return load_gfc(path, nmax);
         },
          py::arg("path"), py::arg("nmax") = -1,
          "Reads an ICGEM .gfc file. nmax < 0 keeps every degree up to max_degree (units: see Model).\n"
          "Raises RuntimeError for a missing or invalid file, ValueError if nmax > max_degree.");
    m.def("load_gfc_sigma", [](const std::string &path, int nmax) {
           py::gil_scoped_release g;
           return load_gfc_sigma(path, nmax);
         },
          py::arg("path"), py::arg("nmax") = -1,
          "Like load_gfc, but C and S hold the formal standard deviations sigma C and sigma S.");

    // ---- Isa and SynthOptions ----------------------------------------------------------------------------------
    py::enum_<Isa>(m, "Isa", "Instruction set for the CPU kernels.")
        .value("Auto", Isa::Auto)
        .value("Scalar", Isa::Scalar)
        .value("Avx2", Isa::Avx2)
        .value("Avx512", Isa::Avx512);

    py::class_<SynthOptions>(m, "SynthOptions", "Synthesis tuning switches; the defaults are the standard path.")
        .def(py::init<>())
        .def_readwrite("threads", &SynthOptions::threads, "OpenMP threads; 0 = OpenMP maximum")
        .def_readwrite("tile_vectors", &SynthOptions::tile_vectors, "SIMD vectors per tile (1, 2, 4, 8, 16); 0 = automatic")
        .def_readwrite("isa", &SynthOptions::isa, "Isa.Auto / Scalar / Avx2 / Avx512")
        .def_readwrite("polar_optimization", &SynthOptions::polar_optimization)
        .def_readwrite("skip_forbidden", &SynthOptions::skip_forbidden)
        .def_readwrite("rotate_trig", &SynthOptions::rotate_trig)
        .def_readwrite("sort_latitude", &SynthOptions::sort_latitude)
        .def_readwrite("parallel_orders", &SynthOptions::parallel_orders, "-1 automatic, 0 per tile, 1 per order m");

    // ---- Points ------------------------------------------------------------------------------------------------
    py::class_<Points>(m, "Points", "Scattered evaluation points. Stored units: lat and lon in RADIANS (geocentric), r in metres.")
        .def(py::init<>())
        .def(py::init([](const DArr &lat, const DArr &lon, const DArr &r) {
               Points p;
               p.lat = copy_in<AlignedVec<double>>(lat, "lat");
               p.lon = copy_in<AlignedVec<double>>(lon, "lon");
               p.r = copy_in<AlignedVec<double>>(r, "r");
               require_same(static_cast<ssize>(p.lat.size()), static_cast<ssize>(p.lon.size()), "lat, lon");
               require_same(static_cast<ssize>(p.lat.size()), static_cast<ssize>(p.r.size()), "lat, r");
               return p;
             }),
             py::arg("lat"), py::arg("lon"), py::arg("r"),
             "Points from geocentric latitude [rad], longitude [rad] and radius [m] (the C++ fields, copied).")
        .def("__len__", &Points::size)
        .def_property_readonly("size", &Points::size)
        .def_property_readonly("lat", [](py::object self) { return view_of(self.cast<const Points &>().lat, self); },
                               "geocentric latitude [rad] (read-only view)")
        .def_property_readonly("lon", [](py::object self) { return view_of(self.cast<const Points &>().lon, self); },
                               "longitude [rad] (read-only view)")
        .def_property_readonly("r", [](py::object self) { return view_of(self.cast<const Points &>().r, self); },
                               "geocentric radius [m] (read-only view)")
        .def_static("from_spherical",
                    [](const DArr &lat_deg, const DArr &lon_deg, const DArr &r) {
                      auto la = copy_in<std::vector<double>>(lat_deg, "lat_deg");
                      auto lo = copy_in<std::vector<double>>(lon_deg, "lon_deg");
                      auto rr = copy_in<std::vector<double>>(r, "r");
                      return Points::from_spherical(la, lo, rr);
                    },
                    py::arg("lat_deg"), py::arg("lon_deg"), py::arg("r"),
                    "Geocentric latitude and longitude in DEGREES, radius in metres.")
        .def_static("from_geodetic",
                    [](const DArr &lat_deg, const DArr &lon_deg, const DArr &h, double a, double e2) {
                      auto la = copy_in<std::vector<double>>(lat_deg, "lat_deg");
                      auto lo = copy_in<std::vector<double>>(lon_deg, "lon_deg");
                      auto hh = copy_in<std::vector<double>>(h, "h");
                      return Points::from_geodetic(la, lo, hh, a, e2);
                    },
                    py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("a") = Ellipsoid::GRS80().a,
                    py::arg("e2") = Ellipsoid::GRS80().e2,
                    "Geodetic latitude and longitude in DEGREES, ellipsoidal height h in metres, on the ellipsoid (a [m], e2).\n"
                    "The stored latitude is converted to geocentric. Defaults: GRS80.");

    // ---- read_table --------------------------------------------------------------------------------------------
    py::class_<PyTable>(m, "Table", "Numeric columns read by read_table.")
        .def_readonly("columns", &PyTable::columns, "list of float64 arrays, one per requested column")
        .def_readonly("skipped_lines", &PyTable::skipped_lines, "1-based numbers of rejected lines (not repaired)");

    m.def("read_table", [](const std::string &path, py::object sep, int header_lines, const std::vector<int> &cols) {
           char delim = '\0'; // '\0' (sep None or empty) = runs of spaces or tabs
           if (!sep.is_none())
           {
             const std::string s = sep.cast<std::string>();
             if (s.size() > 1)
               throw py::value_error("sep must be a single character, or None for whitespace");
             if (s.size() == 1)
               delim = s[0];
           }
           Table t;
           {
             py::gil_scoped_release g;
             t = gsh::read_table(path, delim, header_lines, cols);
           }
           PyTable out;
           for (auto &col : t.columns)
             out.columns.append(to_numpy(std::move(col)));
           out.skipped_lines = std::move(t.skipped_lines);
           return out;
         },
          py::arg("path"), py::arg("sep") = py::none(), py::arg("header_lines") = 0, py::arg("cols"),
          "Reads numeric columns (zero-based indices cols, in output order). sep: one character, or None for whitespace.\n"
          "Lines with a non-numeric requested field are listed in skipped_lines and never repaired.");
  }
}
