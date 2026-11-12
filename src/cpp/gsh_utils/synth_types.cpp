// synth_types.cpp: constructors of Points (moved from gsh_expansion/shs.cpp). They convert latitude
// and longitude from degrees to radians and, for geodetic input, convert to geocentric spherical
// coordinates.
#include "synth_types.hpp"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace gsh
{
  // Degrees to radians for lat and lon; the radius is copied unchanged. Input sizes are checked.
  Points Points::from_spherical(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg, const std::vector<double> &r)
  {
    if (lat_deg.size() != lon_deg.size() || lat_deg.size() != r.size())
      throw std::invalid_argument("Points: input vectors have different sizes");
    Points p;
    p.lat.resize(lat_deg.size());
    p.lon.resize(lat_deg.size());
    p.r.resize(lat_deg.size());
    for (std::size_t i = 0; i < lat_deg.size(); ++i)
    {
      p.lat[i] = deg2rad(lat_deg[i]);
      p.lon[i] = deg2rad(lon_deg[i]);
      p.r[i] = r[i];
    }
    return p;
  }

  // Geodetic input: each point is converted with geodetic_to_spherical (closed form). The stored
  // latitude is geocentric, in radians; the longitude is converted to radians; r is the geocentric
  // radius in metres.
  Points Points::from_geodetic(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg, const std::vector<double> &h, double a, double e2)
  {
    if (lat_deg.size() != lon_deg.size() || lat_deg.size() != h.size())
      throw std::invalid_argument("Points: input vectors have different sizes");
    Points p;
    p.lat.resize(lat_deg.size());
    p.lon.resize(lat_deg.size());
    p.r.resize(lat_deg.size());
    for (std::size_t i = 0; i < lat_deg.size(); ++i)
    {
      const auto s = geodetic_to_spherical(lat_deg[i], lon_deg[i], h[i], a, e2);
      p.lat[i] = s.lat;
      p.lon[i] = s.lon;
      p.r[i] = s.r;
    }
    return p;
  }
}
