// geo.cpp: Ellipsoid (ellipsoid.hpp), coordinate conversions (submodule coords, coordinates.hpp) and permanent-tide
// conversions (submodule tides, tides.hpp, tide_system.hpp).
//
// Array overloads take float64 arrays and call the C++ vector overloads (OpenMP); scalar overloads call the scalar
// C++ functions. Overload selection: a Python float or a NumPy scalar selects the scalar version; a sequence or a
// NumPy array selects the array version.
#include "common.hpp"
#include "coordinates.hpp"
#include "ellipsoid.hpp"
#include "model.hpp"
#include "tide_system.hpp"
#include "tides.hpp"

namespace gshpy
{
  void bind_ellipsoid(py::module_ &m)
  {
    using namespace gsh;

    py::class_<Ellipsoid::Normal>(m, "NormalGravity", "Normal gravity [m/s^2]: g_x north, g_z up, gamma magnitude.")
        .def_readonly("g_x", &Ellipsoid::Normal::g_x)
        .def_readonly("g_z", &Ellipsoid::Normal::g_z)
        .def_readonly("gamma", &Ellipsoid::Normal::gamma);

    py::class_<Ellipsoid>(m, "Ellipsoid",
                          "Reference ellipsoid of revolution. Units: a [m], GM [m^3/s^2], omega [rad/s]; angles in radians\n"
                          "unless the parameter name says _deg. Default: GRS80.")
        .def(py::init<>())
        .def_readwrite("a", &Ellipsoid::a, "semi-major axis [m]")
        .def_readwrite("e2", &Ellipsoid::e2, "first eccentricity squared [-]")
        .def_readwrite("GM", &Ellipsoid::GM, "gravitational constant [m^3/s^2]")
        .def_readwrite("J2", &Ellipsoid::J2, "dynamical form factor [-]")
        .def_readwrite("omega", &Ellipsoid::omega, "angular velocity [rad/s]")
        .def_static("GRS80", &Ellipsoid::GRS80, "the GRS80 ellipsoid")
        .def("zonal", &Ellipsoid::zonal, py::arg("k"), "fully normalised C-bar_{2k,0} of the ellipsoid, k = 0..10")
        .def("radius", &Ellipsoid::radius, py::arg("lat"), "geocentric radius [m] at geocentric latitude lat [rad]")
        .def("geodetic_of_spherical", &Ellipsoid::geodetic_of_spherical, py::arg("lat"),
             "geodetic latitude [rad] from geocentric latitude lat [rad]")
        .def("normal_gravity", &Ellipsoid::normal_gravity, py::arg("lat"), py::arg("r"),
             "normal gravity [m/s^2] at geocentric latitude lat [rad] and radius r [m] (spherical-harmonic series)")
        .def("somigliana", &Ellipsoid::somigliana, py::arg("lat_geodetic"), py::arg("gamma_e") = 9.7803267715,
             py::arg("gamma_p") = 9.8321863685,
             "Somigliana closed-form normal gravity [m/s^2] at geodetic latitude lat_geodetic [rad] (validation)");

    m.def("disturbing_model", &disturbing_model, py::arg("model"), py::arg("ell") = Ellipsoid::GRS80(),
          "Copy of model with the ellipsoidal zonal part removed (delta C-bar_{2k,0}); S unchanged.");
  }

  void bind_coords(py::module_ &m)
  {
    using namespace gsh;
    namespace c = gsh::coords;

    py::class_<c::Ecef>(m, "Ecef", "ECEF position [m].")
        .def(py::init([](double x, double y, double z) { return c::Ecef{x, y, z}; }), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &c::Ecef::x)
        .def_readwrite("y", &c::Ecef::y)
        .def_readwrite("z", &c::Ecef::z);
    py::class_<c::Geodetic>(m, "Geodetic", "Geodetic position: lat_deg, lon_deg [degrees], h [m].")
        .def(py::init([](double lat, double lon, double h) { return c::Geodetic{lat, lon, h}; }), py::arg("lat_deg"),
             py::arg("lon_deg"), py::arg("h"))
        .def_readwrite("lat_deg", &c::Geodetic::lat_deg)
        .def_readwrite("lon_deg", &c::Geodetic::lon_deg)
        .def_readwrite("h", &c::Geodetic::h);
    py::class_<c::SphericalDeg>(m, "SphericalDeg", "Geocentric spherical position: lat_deg, lon_deg [degrees], r [m].")
        .def(py::init([](double lat, double lon, double r) { return c::SphericalDeg{lat, lon, r}; }), py::arg("lat_deg"),
             py::arg("lon_deg"), py::arg("r"))
        .def_readwrite("lat_deg", &c::SphericalDeg::lat_deg)
        .def_readwrite("lon_deg", &c::SphericalDeg::lon_deg)
        .def_readwrite("r", &c::SphericalDeg::r);

    // geodetic_to_ecef: scalar -> Ecef; arrays -> tuple (x, y, z) of arrays.
    m.def("geodetic_to_ecef", [](double lat, double lon, double h, const Ellipsoid &ell) {
           return c::geodetic_to_ecef(lat, lon, h, ell);
         },
          py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("ell") = Ellipsoid::GRS80(),
          "Geodetic (degrees, degrees, metres) to ECEF [m]. Scalars give Ecef; arrays give a tuple (x, y, z).");
    m.def("geodetic_to_ecef", [](const DArr &lat, const DArr &lon, const DArr &h, const Ellipsoid &ell) {
           auto la = copy_in<std::vector<double>>(lat, "lat_deg");
           auto lo = copy_in<std::vector<double>>(lon, "lon_deg");
           auto hh = copy_in<std::vector<double>>(h, "h");
           require_same(static_cast<ssize>(la.size()), static_cast<ssize>(lo.size()), "lat_deg, lon_deg");
           require_same(static_cast<ssize>(la.size()), static_cast<ssize>(hh.size()), "lat_deg, h");
           std::vector<c::Ecef> p;
           {
             py::gil_scoped_release g;
             p = c::geodetic_to_ecef(la, lo, hh, ell);
           }
           std::vector<double> x, y, z;
           for (const auto &q : p)
           {
             x.push_back(q.x);
             y.push_back(q.y);
             z.push_back(q.z);
           }
           return py::make_tuple(to_numpy(std::move(x)), to_numpy(std::move(y)), to_numpy(std::move(z)));
         },
          py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("ell") = Ellipsoid::GRS80());

    // ecef_to_geodetic: Ecef -> Geodetic; arrays x, y, z -> tuple (lat_deg, lon_deg, h).
    m.def("ecef_to_geodetic", [](const c::Ecef &p, const Ellipsoid &ell) { return c::ecef_to_geodetic(p, ell); },
          py::arg("p"), py::arg("ell") = Ellipsoid::GRS80(), "ECEF [m] to geodetic. Returns Geodetic.");
    m.def("ecef_to_geodetic", [](const DArr &x, const DArr &y, const DArr &z, const Ellipsoid &ell) {
           auto xx = copy_in<std::vector<double>>(x, "x");
           auto yy = copy_in<std::vector<double>>(y, "y");
           auto zz = copy_in<std::vector<double>>(z, "z");
           require_same(static_cast<ssize>(xx.size()), static_cast<ssize>(yy.size()), "x, y");
           require_same(static_cast<ssize>(xx.size()), static_cast<ssize>(zz.size()), "x, z");
           std::vector<c::Ecef> p(xx.size());
           for (std::size_t i = 0; i < p.size(); ++i)
             p[i] = c::Ecef{xx[i], yy[i], zz[i]};
           std::vector<c::Geodetic> g;
           {
             py::gil_scoped_release gil;
             g = c::ecef_to_geodetic(p, ell);
           }
           std::vector<double> la, lo, hh;
           for (const auto &q : g)
           {
             la.push_back(q.lat_deg);
             lo.push_back(q.lon_deg);
             hh.push_back(q.h);
           }
           return py::make_tuple(to_numpy(std::move(la)), to_numpy(std::move(lo)), to_numpy(std::move(hh)));
         },
          py::arg("x"), py::arg("y"), py::arg("z"), py::arg("ell") = Ellipsoid::GRS80());

    m.def("geodetic_to_geocentric_latitude", [](double lat, const Ellipsoid &ell) {
           return c::geodetic_to_geocentric_latitude(lat, ell);
         },
          py::arg("lat_deg"), py::arg("ell") = Ellipsoid::GRS80(), "geodetic latitude [deg] -> geocentric latitude [deg]");
    m.def("geodetic_to_geocentric_latitude", [](const DArr &lat, const Ellipsoid &ell) {
           return to_numpy(c::geodetic_to_geocentric_latitude(copy_in<std::vector<double>>(lat, "lat_deg"), ell));
         },
          py::arg("lat_deg"), py::arg("ell") = Ellipsoid::GRS80());

    m.def("geocentric_to_geodetic_latitude", [](double lat, const Ellipsoid &ell) {
           return c::geocentric_to_geodetic_latitude(lat, ell);
         },
          py::arg("lat_deg"), py::arg("ell") = Ellipsoid::GRS80(), "geocentric latitude [deg] -> geodetic latitude [deg]");
    m.def("geocentric_to_geodetic_latitude", [](const DArr &lat, const Ellipsoid &ell) {
           return to_numpy(c::geocentric_to_geodetic_latitude(copy_in<std::vector<double>>(lat, "lat_deg"), ell));
         },
          py::arg("lat_deg"), py::arg("ell") = Ellipsoid::GRS80());

    // geodetic_to_spherical: scalar -> SphericalDeg; arrays -> tuple (lat_deg, lon_deg, r).
    m.def("geodetic_to_spherical", [](double lat, double lon, double h, const Ellipsoid &ell) {
           return c::geodetic_to_spherical(lat, lon, h, ell);
         },
          py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("ell") = Ellipsoid::GRS80(),
          "Geodetic to geocentric spherical (degrees, degrees, metres). Scalars give SphericalDeg; arrays a tuple.");
    m.def("geodetic_to_spherical", [](const DArr &lat, const DArr &lon, const DArr &h, const Ellipsoid &ell) {
           auto la = copy_in<std::vector<double>>(lat, "lat_deg");
           auto lo = copy_in<std::vector<double>>(lon, "lon_deg");
           auto hh = copy_in<std::vector<double>>(h, "h");
           require_same(static_cast<ssize>(la.size()), static_cast<ssize>(lo.size()), "lat_deg, lon_deg");
           require_same(static_cast<ssize>(la.size()), static_cast<ssize>(hh.size()), "lat_deg, h");
           auto s = c::geodetic_to_spherical(la, lo, hh, ell);
           std::vector<double> sla, slo, sr;
           for (const auto &q : s)
           {
             sla.push_back(q.lat_deg);
             slo.push_back(q.lon_deg);
             sr.push_back(q.r);
           }
           return py::make_tuple(to_numpy(std::move(sla)), to_numpy(std::move(slo)), to_numpy(std::move(sr)));
         },
          py::arg("lat_deg"), py::arg("lon_deg"), py::arg("h"), py::arg("ell") = Ellipsoid::GRS80());
  }

  void bind_tides(py::module_ &m)
  {
    using namespace gsh;
    namespace t = gsh::tides;

    // ---- constants --------------------------------------------------------------------------------------------
    m.attr("A0") = t::A0;
    m.attr("H0") = t::H0;
    m.attr("GM") = t::GM;
    m.attr("R") = t::R;
    m.attr("GAMMA_E") = t::GAMMA_E;
    m.attr("GRS80_K") = t::GRS80_K;
    m.attr("GRS80_E2") = t::GRS80_E2;
    m.attr("K20_DEFAULT") = t::K20_DEFAULT;
    m.attr("NO_C20_MT") = std::string(t::NO_C20_MT);

    // ---- tide systems -----------------------------------------------------------------------------------------
    py::enum_<TideSystem>(m, "TideSystem", "Permanent-tide system: TideFree (NT), ZeroTide (ZT), MeanTide (MT).")
        .value("TideFree", TideSystem::TideFree)
        .value("ZeroTide", TideSystem::ZeroTide)
        .value("MeanTide", TideSystem::MeanTide);
    m.def("parse_tide", [](const std::string &name) { return parse_tide(name); }, py::arg("name"),
          "Parses a system name ('tide_free', 'zero_tide', 'mean_tide' and aliases). Raises ValueError otherwise.");
    m.def("to_string", [](TideSystem v) { return std::string(to_string(v)); }, py::arg("system"),
          "Canonical name of a TideSystem.");

    // ---- C20 and shifts ---------------------------------------------------------------------------------------
    m.def("c20_shift", [](TideSystem from, TideSystem to, double k2) { return t::c20_shift(from, to, k2); },
          py::arg("from_system"), py::arg("to_system"), py::arg("k2") = t::K20_DEFAULT,
          "Change of C-bar20 (to minus from). Raises ValueError for MeanTide (no C20 in mean-tide).");

    // height_shift: N shift [m] at lat (degrees).
    m.def("height_shift", [](double lat, TideSystem from, TideSystem to, double k) { return t::height_shift(lat, from, to, k); },
          py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), py::arg("k") = t::K20_DEFAULT,
          "Amount [m] to ADD to the geoid height N, at geodetic latitude lat_deg [degrees].");
    m.def("height_shift", [](const DArr &lat, TideSystem from, TideSystem to, double k) {
           return to_numpy(t::height_shift(copy_in<std::vector<double>>(lat, "lat_deg"), from, to, k));
         },
          py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), py::arg("k") = t::K20_DEFAULT);

    m.def("convert_geoid_height", [](double N, double lat, TideSystem from, TideSystem to, double k) {
           return t::convert_geoid_height(N, lat, from, to, k);
         },
          py::arg("N"), py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), py::arg("k") = t::K20_DEFAULT,
          "Geoid height N [m] converted between systems.");
    m.def("convert_geoid_height", [](const DArr &N, const DArr &lat, TideSystem from, TideSystem to, double k) {
           auto nn = copy_in<std::vector<double>>(N, "N");
           auto la = copy_in<std::vector<double>>(lat, "lat_deg");
           require_same(static_cast<ssize>(nn.size()), static_cast<ssize>(la.size()), "N, lat_deg");
           return to_numpy(t::convert_geoid_height(nn, la, from, to, k));
         },
          py::arg("N"), py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), py::arg("k") = t::K20_DEFAULT);

    m.def("convert_model_c20", [](const Model &model, TideSystem to, double k2) { return t::convert_model_c20(model, to, k2); },
          py::arg("model"), py::arg("to_system"), py::arg("k2") = t::K20_DEFAULT,
          "Copy of model with C-bar20 moved to system to_system (source = model.tide_system).");
    m.def("convert_c20_from", [](const Model &model, TideSystem from, TideSystem to, double k2) {
           return t::convert_c20_from(model, from, to, k2);
         },
          py::arg("model"), py::arg("from_system"), py::arg("to_system"), py::arg("k2") = t::K20_DEFAULT,
          "Copy of model with C-bar20 converted from an explicit source system.");

    // ---- lat-only terms (degrees in, metres or mGal out) ----------------------------------------------------
    m.def("crust_height_shift_m", [](double lat) { return t::crust_height_shift_m(lat); }, py::arg("lat_deg"),
          "Vertical crust displacement h_T [m] at lat_deg [degrees].");
    m.def("crust_height_shift_m", [](const DArr &lat) {
           return to_numpy(t::crust_height_shift_m(copy_in<std::vector<double>>(lat, "lat_deg")));
         },
          py::arg("lat_deg"));
    m.def("geoid_tide_separation_m", [](double lat) { return t::geoid_tide_separation_m(lat); }, py::arg("lat_deg"),
          "Mean-minus-zero geoid separation H_T [m] at lat_deg [degrees].");
    m.def("geoid_tide_separation_m", [](const DArr &lat) {
           return to_numpy(t::geoid_tide_separation_m(copy_in<std::vector<double>>(lat, "lat_deg")));
         },
          py::arg("lat_deg"));
    m.def("gravity_tide_mgal", [](double lat) { return t::gravity_tide_mgal(lat); }, py::arg("lat_deg"),
          "Direct tidal gravity term g_T [mGal] at lat_deg [degrees].");
    m.def("gravity_tide_mgal", [](const DArr &lat) {
           return to_numpy(t::gravity_tide_mgal(copy_in<std::vector<double>>(lat, "lat_deg")));
         },
          py::arg("lat_deg"));
    m.def("gravity_tide_shift_mgal", [](double lat, TideSystem from, TideSystem to) {
           return t::gravity_tide_shift_mgal(lat, from, to);
         },
          py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), "Amount [mGal] to ADD to observed gravity.");
    m.def("gravity_tide_shift_mgal", [](const DArr &lat, TideSystem from, TideSystem to) {
           return to_numpy(t::gravity_tide_shift_mgal(copy_in<std::vector<double>>(lat, "lat_deg"), from, to));
         },
          py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"));

    m.def("convert_ellipsoidal_height", [](double h, double lat, TideSystem from, TideSystem to) {
           return t::convert_ellipsoidal_height(h, lat, from, to);
         },
          py::arg("h"), py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"), "Ellipsoidal height h [m] converted.");
    m.def("convert_ellipsoidal_height", [](const DArr &h, const DArr &lat, TideSystem from, TideSystem to) {
           auto hh = copy_in<std::vector<double>>(h, "h");
           auto la = copy_in<std::vector<double>>(lat, "lat_deg");
           require_same(static_cast<ssize>(hh.size()), static_cast<ssize>(la.size()), "h, lat_deg");
           return to_numpy(t::convert_ellipsoidal_height(hh, la, from, to));
         },
          py::arg("h"), py::arg("lat_deg"), py::arg("from_system"), py::arg("to_system"));

    m.def("ihrf_separation_surface", [](double zeta, double lat) { return t::ihrf_separation_surface(zeta, lat); },
          py::arg("zeta_zero_tide"), py::arg("lat_deg"), "Separation surface S [m] (IHRF).");
    m.def("ihrf_separation_surface", [](const DArr &zeta, const DArr &lat) {
           auto zz = copy_in<std::vector<double>>(zeta, "zeta_zero_tide");
           auto la = copy_in<std::vector<double>>(lat, "lat_deg");
           require_same(static_cast<ssize>(zz.size()), static_cast<ssize>(la.size()), "zeta_zero_tide, lat_deg");
           return to_numpy(t::ihrf_separation_surface(zz, la));
         },
          py::arg("zeta_zero_tide"), py::arg("lat_deg"));
  }
}
