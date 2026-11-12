// test_coordinates.cpp: tests of gsh::coords (geodetic <-> ECEF <-> geocentric conversions) on GRS80.
// Covers exact closed-form values (equator, poles), known numbers (geocentric latitude at 45 deg),
// round trips over a grid that includes the poles, the equator, h = 0, h = 8848 m and h = 1e5 m,
// consistency with gsh::geodetic_to_spherical, and that the vector overloads reproduce the scalar
// functions exactly.
#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "coordinates.hpp"
#include "core.hpp"
#include "ellipsoid.hpp"

namespace
{
  using gsh::Ellipsoid;
  namespace coords = gsh::coords;

  const Ellipsoid G = Ellipsoid::GRS80();
  // Semi-minor axis of GRS80: b = a sqrt(1 - e2) = 6356752.31414036 m.
  const double B = G.a * std::sqrt(1.0 - G.e2);

  // Smallest angular difference between two longitudes, in degrees (handles the +-180 wrap).
  double lon_diff_deg(double a, double b)
  {
    return std::fabs(std::remainder(a - b, 360.0));
  }
} // namespace

// Equator, lambda = 0, h = 0: the point lies on the x axis at distance a.
TEST(Coordinates, EquatorAtZeroLongitudeIsSemiMajorAxis)
{
  const auto e = coords::geodetic_to_ecef(0.0, 0.0, 0.0, G);
  EXPECT_NEAR(e.x, G.a, 1e-9);
  EXPECT_NEAR(e.y, 0.0, 1e-9);
  EXPECT_NEAR(e.z, 0.0, 1e-9);
}

// Equator, lambda = 90 deg: the point lies on the y axis.
TEST(Coordinates, EquatorAtNinetyEastIsOnYAxis)
{
  const auto e = coords::geodetic_to_ecef(0.0, 90.0, 0.0, G);
  EXPECT_NEAR(e.x, 0.0, 1e-9);
  EXPECT_NEAR(e.y, G.a, 1e-9);
  EXPECT_NEAR(e.z, 0.0, 1e-9);
}

// North pole, h = 0: the point lies on the z axis at the semi-minor axis b.
TEST(Coordinates, NorthPoleIsSemiMinorAxis)
{
  const auto e = coords::geodetic_to_ecef(90.0, 0.0, 0.0, G);
  EXPECT_NEAR(e.x, 0.0, 1e-9);
  EXPECT_NEAR(e.y, 0.0, 1e-9);
  EXPECT_NEAR(e.z, B, 1e-6);
  EXPECT_NEAR(e.z, 6356752.3141, 1e-3); // published GRS80 polar radius
}

// South pole, h = 0: symmetric to the north pole.
TEST(Coordinates, SouthPoleIsMinusSemiMinorAxis)
{
  const auto e = coords::geodetic_to_ecef(-90.0, 0.0, 0.0, G);
  EXPECT_NEAR(e.z, -B, 1e-6);
}

// Equator with height: X = a + h exactly (the normal is the x axis).
TEST(Coordinates, EquatorHeightAddsToRadius)
{
  const auto e = coords::geodetic_to_ecef(0.0, 0.0, 8848.0, G);
  EXPECT_NEAR(e.x, G.a + 8848.0, 1e-9);
  EXPECT_NEAR(e.z, 0.0, 1e-9);
}

// Inverse at the exact pole and equator (p = 0 and z = 0 are the singular cases of the iteration).
TEST(Coordinates, InverseAtPolesAndEquator)
{
  const auto north = coords::ecef_to_geodetic(coords::Ecef{0.0, 0.0, B + 100.0}, G);
  EXPECT_NEAR(north.lat_deg, 90.0, 1e-12);
  EXPECT_NEAR(north.h, 100.0, 1e-7);

  const auto south = coords::ecef_to_geodetic(coords::Ecef{0.0, 0.0, -B - 50.0}, G);
  EXPECT_NEAR(south.lat_deg, -90.0, 1e-12);
  EXPECT_NEAR(south.h, 50.0, 1e-7);

  const auto eq = coords::ecef_to_geodetic(coords::Ecef{G.a + 8848.0, 0.0, 0.0}, G);
  EXPECT_NEAR(eq.lat_deg, 0.0, 1e-12);
  EXPECT_NEAR(eq.lon_deg, 0.0, 1e-12);
  EXPECT_NEAR(eq.h, 8848.0, 1e-9);
}

// Geocentric latitude at 45 deg geodetic: tan(phi') = (1 - e2) tan(45 deg). The reference value
// 44.8075767830733 deg was computed independently as atan((1 - e2) * 1) in double precision.
TEST(Coordinates, GeocentricLatitudeKnownValue)
{
  EXPECT_NEAR(coords::geodetic_to_geocentric_latitude(45.0, G), 44.8075767830733, 1e-10);
  EXPECT_NEAR(coords::geocentric_to_geodetic_latitude(44.8075767830733, G), 45.0, 1e-10);
}

// The geocentric latitude equals the geodetic one at the equator and at the poles.
TEST(Coordinates, GeocentricLatitudeFixedPoints)
{
  EXPECT_DOUBLE_EQ(coords::geodetic_to_geocentric_latitude(0.0, G), 0.0);
  EXPECT_NEAR(coords::geodetic_to_geocentric_latitude(90.0, G), 90.0, 1e-12);
  EXPECT_NEAR(coords::geodetic_to_geocentric_latitude(-90.0, G), -90.0, 1e-12);
  EXPECT_NEAR(coords::geocentric_to_geodetic_latitude(90.0, G), 90.0, 1e-12);
}

// Round trip geodetic -> ECEF -> geodetic on a grid that includes the poles, the equator, h = 0,
// h = 8848 m and h = 1e5 m. Accuracy observed with this implementation: |dh| below 1e-9 m, latitude
// and longitude below 1e-13 deg. Longitude is not compared at the poles, where it is undefined.
TEST(Coordinates, RoundTripGrid)
{
  const double lats[] = {-90.0, -89.999999, -60.0, -45.0, -1e-6, 0.0, 1e-6, 0.5, 30.0, 45.0, 89.999999, 90.0};
  const double lons[] = {-180.0, -90.0, -1e-3, 0.0, 45.0, 90.0, 179.9};
  const double hs[] = {-1e4, -1.0, 0.0, 8848.0, 1e5};
  for (double la : lats)
    for (double lo : lons)
      for (double h : hs)
      {
        const auto e = coords::geodetic_to_ecef(la, lo, h, G);
        const auto d = coords::ecef_to_geodetic(e, G);
        const std::string tag = "lat=" + std::to_string(la) + " lon=" + std::to_string(lo) + " h=" + std::to_string(h);
        EXPECT_NEAR(d.h, h, 1e-9) << tag;
        EXPECT_NEAR(d.lat_deg, la, 1e-11) << tag;
        if (std::fabs(la) < 89.9999)
          EXPECT_LT(lon_diff_deg(d.lon_deg, lo), 1e-11) << tag;
      }
}

// Round trip for the summit of Mount Everest and for a point at h = 0 in the southern hemisphere.
TEST(Coordinates, RoundTripKnownPoints)
{
  const auto everest = coords::ecef_to_geodetic(coords::geodetic_to_ecef(27.9881, 86.9250, 8848.0, G), G);
  EXPECT_NEAR(everest.lat_deg, 27.9881, 1e-11);
  EXPECT_NEAR(everest.lon_deg, 86.9250, 1e-11);
  EXPECT_NEAR(everest.h, 8848.0, 1e-9);

  const auto south = coords::ecef_to_geodetic(coords::geodetic_to_ecef(-33.45, -70.66, 0.0, G), G);
  EXPECT_NEAR(south.lat_deg, -33.45, 1e-11);
  EXPECT_NEAR(south.lon_deg, -70.66, 1e-11);
  EXPECT_NEAR(south.h, 0.0, 1e-9);
}

// coords::geodetic_to_spherical must agree with gsh::geodetic_to_spherical (core.hpp), which is the
// shared implementation. Units differ: the core version returns radians.
TEST(Coordinates, SphericalAgreesWithCore)
{
  const double lats[] = {-80.0, -10.0, 0.0, 45.0, 89.0};
  const double lons[] = {-120.0, 0.0, 77.0};
  for (double la : lats)
    for (double lo : lons)
    {
      const auto s = coords::geodetic_to_spherical(la, lo, 8848.0, G);
      const auto c = gsh::geodetic_to_spherical(la, lo, 8848.0, G.a, G.e2);
      EXPECT_NEAR(gsh::deg2rad(s.lat_deg), c.lat, 1e-15);
      EXPECT_DOUBLE_EQ(s.lon_deg, lo);
      EXPECT_DOUBLE_EQ(s.r, c.r);
    }
}

// Spherical radius at the poles and at the equator (h = 0).
TEST(Coordinates, SphericalRadiusAtPolesAndEquator)
{
  EXPECT_NEAR(coords::geodetic_to_spherical(90.0, 0.0, 0.0, G).r, B, 1e-6);
  EXPECT_NEAR(coords::geodetic_to_spherical(0.0, 0.0, 0.0, G).r, G.a, 1e-9);
  EXPECT_NEAR(coords::geodetic_to_spherical(90.0, 0.0, 0.0, G).lat_deg, 90.0, 1e-12);
}

// The vector overloads run the scalar code element by element, so the results must be identical.
TEST(Coordinates, VectorOverloadsEqualScalar)
{
  const std::vector<double> lat = {-90.0, -45.0, -1e-3, 0.0, 12.5, 45.0, 89.9, 90.0};
  const std::vector<double> lon = {-179.0, -90.0, 0.0, 33.3, 90.0, 120.0, 180.0, 0.0};
  const std::vector<double> h = {-1e4, 0.0, 8848.0, 1e5, 10.0, -5.0, 0.0, 250.0};

  const auto ecef = coords::geodetic_to_ecef(lat, lon, h, G);
  const auto geo = coords::ecef_to_geodetic(ecef, G);
  const auto gc = coords::geodetic_to_geocentric_latitude(lat, G);
  const auto gd = coords::geocentric_to_geodetic_latitude(gc, G);
  const auto sph = coords::geodetic_to_spherical(lat, lon, h, G);
  ASSERT_EQ(ecef.size(), lat.size());
  ASSERT_EQ(geo.size(), lat.size());
  ASSERT_EQ(gc.size(), lat.size());
  ASSERT_EQ(gd.size(), lat.size());
  ASSERT_EQ(sph.size(), lat.size());

  for (std::size_t i = 0; i < lat.size(); ++i)
  {
    const auto e = coords::geodetic_to_ecef(lat[i], lon[i], h[i], G);
    EXPECT_EQ(ecef[i].x, e.x) << i;
    EXPECT_EQ(ecef[i].y, e.y) << i;
    EXPECT_EQ(ecef[i].z, e.z) << i;

    const auto g = coords::ecef_to_geodetic(e, G);
    EXPECT_EQ(geo[i].lat_deg, g.lat_deg) << i;
    EXPECT_EQ(geo[i].lon_deg, g.lon_deg) << i;
    EXPECT_EQ(geo[i].h, g.h) << i;

    EXPECT_EQ(gc[i], coords::geodetic_to_geocentric_latitude(lat[i], G)) << i;
    EXPECT_EQ(gd[i], coords::geocentric_to_geodetic_latitude(gc[i], G)) << i;

    const auto s = coords::geodetic_to_spherical(lat[i], lon[i], h[i], G);
    EXPECT_EQ(sph[i].lat_deg, s.lat_deg) << i;
    EXPECT_EQ(sph[i].lon_deg, s.lon_deg) << i;
    EXPECT_EQ(sph[i].r, s.r) << i;
  }

  // Vector round trip on the same data.
  for (std::size_t i = 0; i < lat.size(); ++i)
  {
    if (std::fabs(lat[i]) < 89.9999)
      EXPECT_LT(lon_diff_deg(geo[i].lon_deg, lon[i]), 1e-11) << i;
    EXPECT_NEAR(geo[i].h, h[i], 1e-9) << i;
  }
}

// Mismatched input vector sizes are rejected; empty input gives empty output.
TEST(Coordinates, VectorOverloadsRejectMismatchedSizes)
{
  EXPECT_THROW(coords::geodetic_to_ecef(std::vector<double>{0.0, 1.0}, std::vector<double>{0.0}, std::vector<double>{0.0, 0.0}, G),
               std::invalid_argument);
  EXPECT_THROW(coords::geodetic_to_spherical(std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{}, G),
               std::invalid_argument);
  EXPECT_TRUE(coords::geodetic_to_ecef(std::vector<double>{}, std::vector<double>{}, std::vector<double>{}, G).empty());
  EXPECT_TRUE(coords::ecef_to_geodetic(std::vector<coords::Ecef>{}, G).empty());
}
