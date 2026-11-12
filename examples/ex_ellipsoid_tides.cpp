// ex_ellipsoid_tides.cpp: GRS80 normal gravity, coordinate conversions and permanent-tide conversions.
//
// Shows: Somigliana's closed form against the spherical-harmonic normal gravity of the ellipsoid (ellipsoid.hpp);
// geodetic, ECEF and geocentric spherical coordinates (coordinates.hpp, core.hpp); and the conversion of C20, of the
// height N, of the ellipsoidal height and of gravity between the tide-free, zero-tide and mean-tide systems (tides.hpp).
//
// Usage: ex_ellipsoid_tides <model.gfc>
//   The model file is read only for its C20 (degree 2). Its header must name the tide system (EGM2008: tide_free).
//
// Output: normal-gravity differences in mGal (1 mGal = 1e-5 m/s^2), round-trip errors of the coordinate conversions, and
// the tide shifts. A shift is defined as "value in the target system minus value in the source system".
// C20 cannot be converted to mean-tide: tides.hpp refuses that request, and the example shows the refusal.

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <stdexcept>
#include "coordinates.hpp"
#include "core.hpp"
#include "ellipsoid.hpp"
#include "model.hpp"
#include "tide_system.hpp"
#include "tides.hpp"

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_ellipsoid_tides: %s\nusage: ex_ellipsoid_tides <model.gfc>\n", why);
    std::exit(2);
  }

  void normal_gravity_section(const gsh::Ellipsoid &ell)
  {
    std::printf("Normal gravity on GRS80, at the surface (h = 0):\n");
    for (double lat : {0.0, 45.0, 80.0})
    {
      const double lat_c = gsh::coords::geodetic_to_geocentric_latitude(lat, ell);
      const double r = ell.radius(gsh::deg2rad(lat_c));
      const double g_sh = ell.normal_gravity(gsh::deg2rad(lat_c), r).gamma;
      const double g_closed = ell.somigliana(gsh::deg2rad(lat));
      std::printf("  lat %5.1f deg: Somigliana %.10f m/s2, spherical harmonics %.10f m/s2, difference %.3e mGal\n", lat,
                  g_closed, g_sh, (g_sh - g_closed) * 1e5);
    }
  }

  void coordinates_section(const gsh::Ellipsoid &ell)
  {
    const double lat = 45.5, lon = -70.25, h = 123.4; // geodetic, degrees and metres
    const gsh::coords::Ecef p = gsh::coords::geodetic_to_ecef(lat, lon, h, ell);
    const gsh::coords::Geodetic back = gsh::coords::ecef_to_geodetic(p, ell);
    const gsh::coords::SphericalDeg s = gsh::coords::geodetic_to_spherical(lat, lon, h, ell);
    const gsh::Spherical core = gsh::geodetic_to_spherical(lat, lon, h, ell.a, ell.e2); // radians
    std::printf("Coordinates of geodetic (%.2f, %.2f, %.1f m):\n", lat, lon, h);
    std::printf("  ECEF: X = %.4f, Y = %.4f, Z = %.4f m\n", p.x, p.y, p.z);
    std::printf("  ECEF back to geodetic: lat %.10f, lon %.10f, h %.6f m (error in h %.2e m)\n", back.lat_deg, back.lon_deg,
                back.h, back.h - h);
    std::printf("  geocentric spherical: lat %.10f deg, lon %.6f deg, r %.4f m\n", s.lat_deg, s.lon_deg, s.r);
    std::printf("  core.hpp version (radians): difference in latitude %.2e rad, in r %.2e m\n",
                core.lat - gsh::deg2rad(s.lat_deg), core.r - s.r);
  }

  void tides_section(const gsh::Model &model)
  {
    using gsh::TideSystem;
    namespace tides = gsh::tides;
    const double lat = 45.0, N = 30.0, h = 100.0; // degrees, metres, metres
    const gsh::Model zt = tides::convert_model_c20(model, TideSystem::ZeroTide);

    std::printf("Permanent-tide systems (file header tide_system = \"%s\", k20 = %.2f):\n", model.tide_system.c_str(),
                tides::K20_DEFAULT);
    std::printf("  C20 tide-free  = %.12e\n", model.c(2, 0));
    std::printf("  C20 zero-tide  = %.12e, shift %.6e\n", zt.c(2, 0), zt.c(2, 0) - model.c(2, 0));
    std::printf("  c20_shift(tide-free -> zero-tide) = %.6e\n", tides::c20_shift(TideSystem::TideFree, TideSystem::ZeroTide));
    try
    {
      tides::c20_shift(TideSystem::TideFree, TideSystem::MeanTide);
    }
    catch (const std::invalid_argument &e)
    {
      std::printf("  c20_shift(tide-free -> mean-tide) refused: %s\n", e.what());
    }

    const TideSystem TF = TideSystem::TideFree, ZT = TideSystem::ZeroTide, MT = TideSystem::MeanTide;
    std::printf("Shifts at lat %.1f deg (value in target minus value in source):\n", lat);
    std::printf("  height N, tide-free -> zero-tide: %.6f m, tide-free -> mean-tide: %.6f m\n",
                tides::height_shift(lat, TF, ZT), tides::height_shift(lat, TF, MT));
    std::printf("  N = %.1f m converted to mean-tide: %.6f m\n", N, tides::convert_geoid_height(N, lat, TF, MT));
    std::printf("  ellipsoidal h = %.1f m converted to mean-tide: %.6f m\n", h,
                tides::convert_ellipsoidal_height(h, lat, TF, MT));
    std::printf("  gravity: direct tidal term g_T = %.6f mGal, tide-free -> mean-tide shift = %.6f mGal\n",
                tides::gravity_tide_mgal(lat), tides::gravity_tide_shift_mgal(lat, TF, MT));
    std::printf("  tide_from_name(\"zero\") = %s\n", gsh::to_string(tides::tide_from_name("zero")));
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    if (argc != 2)
      usage("expected one argument: the model file");
    const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();
    const gsh::Model model = gsh::load_gfc(argv[1], 2);
    normal_gravity_section(ell);
    coordinates_section(ell);
    tides_section(model);
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "ex_ellipsoid_tides: %s\n", e.what());
    return 1;
  }
}
