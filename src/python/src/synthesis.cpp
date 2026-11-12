// synthesis.cpp: raw-sum synthesis at points (shs.hpp) and on latitude-longitude grids (grid.hpp).
//
// Raw sums are dimensionless (see shs.hpp). Field order (RawFields.f) is given by RAW_FIELD_NAMES.
#include <memory>
#include "common.hpp"
#include "grid.hpp"
#include "legendre.hpp"
#include "model.hpp"
#include "shs.hpp"
#include "synth_types.hpp"

namespace gshpy
{
  namespace
  {
    // Python-facing holder for the shared recurrence tables (shared_ptr<const Recurrence> is not a valid holder).
    struct PyRecurrence
    {
      std::shared_ptr<const gsh::Recurrence> rec;
    };
  }

  void bind_synthesis(py::module_ &m)
  {
    using namespace gsh;

    m.attr("RAW_FIELD_NAMES") = py::make_tuple("S0", "SR", "SP", "SL", "SRR", "SRP", "SRL", "SPP", "SPL", "SLL");

    // ---- recurrence and preparation ---------------------------------------------------------------------------
    py::class_<PyRecurrence>(m, "Recurrence", "Recurrence tables for one maximum degree; shareable between Prepared objects.")
        .def_property_readonly("nmax", [](const PyRecurrence &r) { return r.rec->nmax; });
    m.def("make_recurrence", [](int nmax) {
           PyRecurrence r;
           {
             py::gil_scoped_release g;
             r.rec = make_recurrence(nmax);
           }
           return r;
         },
          py::arg("nmax"));

    py::class_<Prepared>(m, "Prepared", "A model prepared for the kernels: coefficients in [nmin, nmax].")
        .def_readonly("nmin", &Prepared::nmin)
        .def_readonly("nmax", &Prepared::nmax)
        .def_readonly("GM", &Prepared::GM)
        .def_readonly("R", &Prepared::R);
    m.def("prepare", [](const Model &model, int nmin, int nmax, py::object recurrence) {
           std::shared_ptr<const Recurrence> rec;
           if (!recurrence.is_none())
             rec = recurrence.cast<const PyRecurrence &>().rec;
           Prepared p;
           {
             py::gil_scoped_release g;
             p = prepare(model, nmin, nmax, rec);
           }
           return p;
         },
          py::arg("model"), py::arg("nmin") = 0, py::arg("nmax") = -1, py::arg("recurrence") = py::none(),
          "Prepares model for synthesis. nmax < 0 means model.nmax. recurrence: a Recurrence for the same nmax, or None.");

    // ---- raw sums ----------------------------------------------------------------------------------------------
    py::class_<RawFields>(m, "RawFields", "Raw sums per point or per grid node; f[k] follows RAW_FIELD_NAMES.")
        .def_readonly("grad", &RawFields::grad)
        .def_property_readonly("f", [](py::object self) {
          const auto &r = self.cast<const RawFields &>();
          py::list out;
          for (const auto &v : r.f)
            out.append(view_of(v, self));
          return out;
        }, "list of read-only float64 arrays, one per field");

    m.def("synth_points", [](const Prepared &p, const Points &pts, int grad, const SynthOptions &o) {
           RawFields r;
           {
             py::gil_scoped_release g;
             r = synth_points(p, pts, grad, o);
           }
           return r;
         },
          py::arg("prep"), py::arg("pts"), py::arg("grad") = 0, py::arg("opts") = SynthOptions(),
          "Raw sums at every point (grad 0, 1 or 2). Dimensionless. Multiply S0 by GM/r to get V.");

    m.def("active_isa_name", [](Isa isa) { return active_isa_name(isa); }, py::arg("isa") = Isa::Auto,
          "Name of the ISA synth_points would use: 'avx512', 'avx2' or 'scalar'.");

    // ---- oracle ------------------------------------------------------------------------------------------------
    py::class_<ReferencePoint>(m, "ReferencePoint", "Full-triangle oracle result at one point.")
        .def_property_readonly("raw", [](const ReferencePoint &r) {
          py::list l;
          for (double x : r.raw)
            l.append(x);
          return l;
        }, "raw sums, RAW_FIELD_NAMES order")
        .def_readonly("vxx", &ReferencePoint::vxx)
        .def_readonly("vxy", &ReferencePoint::vxy)
        .def_readonly("vxz", &ReferencePoint::vxz)
        .def_readonly("vyy", &ReferencePoint::vyy)
        .def_readonly("vyz", &ReferencePoint::vyz)
        .def_readonly("vzz", &ReferencePoint::vzz);

    m.def("reference_point", [](const Model &model, const legendre::Coefs &coefs, double lat, double lon, double r, int nmin,
                                int nmax) {
           const int eff = nmax < 0 ? model.nmax : nmax;
           if (eff > model.nmax)
             throw py::value_error("nmax is larger than the model degree");
           if (coefs.nmax < eff)
             throw py::value_error("coefs.nmax must be at least nmax (use make_coefs(nmax))");
           ReferencePoint rp;
           {
             py::gil_scoped_release g;
             rp = reference_point(model, coefs, lat, lon, r, nmin, nmax);
           }
           return rp;
         },
          py::arg("model"), py::arg("coefs"), py::arg("lat"), py::arg("lon"), py::arg("r"), py::arg("nmin") = 0,
          py::arg("nmax") = -1,
          "Oracle at one point. lat and lon in RADIANS (geocentric), r in metres. Slow: O(nmax^2).");

    // ---- grids -------------------------------------------------------------------------------------------------
    py::class_<Grid>(m, "Grid", "Latitude-longitude grid. lat [rad, north to south], r [m] per row, lon0 and dlon [rad].")
        .def(py::init([](const DArr &lat, const DArr &r, double lon0, double dlon, ssize nlon) {
               Grid g;
               g.lat = copy_in<AlignedVec<double>>(lat, "lat");
               g.r = copy_in<AlignedVec<double>>(r, "r");
               require_same(static_cast<ssize>(g.lat.size()), static_cast<ssize>(g.r.size()), "lat, r");
               g.lon0 = lon0;
               g.dlon = dlon;
               g.nlon = static_cast<std::int64_t>(nlon);
               return g;
             }),
             py::arg("lat"), py::arg("r"), py::arg("lon0"), py::arg("dlon"), py::arg("nlon"),
             "Custom grid: lat [rad], r [m] per row, first longitude lon0 [rad], step dlon [rad], nlon columns.")
        .def_property_readonly("lat", [](py::object self) { return view_of(self.cast<const Grid &>().lat, self); })
        .def_property_readonly("r", [](py::object self) { return view_of(self.cast<const Grid &>().r, self); })
        .def_readonly("lon0", &Grid::lon0)
        .def_readonly("dlon", &Grid::dlon)
        .def_readonly("nlon", &Grid::nlon)
        .def_property_readonly("nlat", [](const Grid &g) { return g.nlat(); })
        .def("full_circle", &Grid::full_circle, "True when nlon * dlon = 2 pi (selects the FFT path).")
        .def_static("gauss_legendre", [](int nmax, double r) {
           py::gil_scoped_release g;
           return Grid::gauss_legendre(nmax, r);
         },
          py::arg("nmax"), py::arg("r"), "Gauss-Legendre grid: nmax+1 latitudes, 2(nmax+1) longitudes. r [m].")
        .def_static("driscoll_healy1", [](int nmax, double r) {
           py::gil_scoped_release g;
           return Grid::driscoll_healy1(nmax, r);
         },
          py::arg("nmax"), py::arg("r"), "Driscoll-Healy grid, variant 1: 2L x 2L.")
        .def_static("driscoll_healy2", [](int nmax, double r) {
           py::gil_scoped_release g;
           return Grid::driscoll_healy2(nmax, r);
         },
          py::arg("nmax"), py::arg("r"), "Driscoll-Healy grid, variant 2: 2L x 4L.")
        .def_static("regular", &Grid::regular, py::arg("lat_north_deg"), py::arg("lat_south_deg"), py::arg("nlat"),
                    py::arg("lon_west_deg"), py::arg("lon_east_deg"), py::arg("nlon"), py::arg("r"),
                    "Regular regional or global grid in DEGREES (both ends included). r [m].");

    py::class_<GridTiming>(m, "GridTiming", "Timing of synth_grid: lumps and fft in thread-seconds, total in wall seconds.")
        .def_readonly("lumps", &GridTiming::lumps)
        .def_readonly("fft", &GridTiming::fft)
        .def_readonly("total", &GridTiming::total);

    m.def("synth_grid", [](const Prepared &p, const Grid &g, int grad, const SynthOptions &o) {
           RawFields r;
           GridTiming t;
           {
             py::gil_scoped_release gil;
             r = synth_grid(p, g, grad, o, &t);
           }
           return py::make_tuple(r, t);
         },
          py::arg("prep"), py::arg("grid"), py::arg("grad") = 0, py::arg("opts") = SynthOptions(),
          "Raw sums on a grid. Returns (RawFields, GridTiming); each field holds nlat*nlon values, latitude slow.");

    m.def("grid_to_points", [](const Grid &g) { return grid_to_points(g); }, py::arg("grid"),
          "Expands a grid into scattered Points (row by row). Slow for large grids.");
  }
}
