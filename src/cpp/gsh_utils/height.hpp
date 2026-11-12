// height.hpp: source of the topographic height H(phi, lambda) used by the N and zeta functionals (functionals.hpp).
// The height enters N through the term -2 pi G rho H^2 / gamma(r_ell, phi) with G = 6.67259e-11 m^3
// kg^-1 s^-2 and rho = 2670 kg m^-3 (Bucha and Janák 2019).
//
// This header provides:
//   - HeightGrid: a regular latitude/longitude grid of heights with bilinear interpolation and block
//     averaging;
//   - load_gravsoft_grid: reader for the GRAVSOFT grid format;
//   - HeightSource: the abstract interface used by the computation, with a grid implementation
//     (height_from_grid).
// The implementation based on a spherical-harmonic topography model (which uses the synthesis code of
// gsh_expansion) is in gsh_expansion/height_model.hpp.
//
// Units: heights in metres, angles in degrees for grid geometry, longitudes in 0..360.
#pragma once
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include "core.hpp"
#include "synth_types.hpp"

namespace gsh
{
  // Regular grid of heights. Rows run from north to south (row 0 = lat_north); columns run from west to
  // east (column 0 = lon_west). Element (y, x) is h[y * nlon + x].
  struct HeightGrid
  {
    double lat_south = 0, lat_north = 0, lon_west = 0, lon_east = 0; // degrees; longitudes 0..360
    int nlat = 0, nlon = 0;   // number of nodes in latitude and longitude
    std::vector<float> h;     // nlat * nlon heights [m]

    // Grid spacing in degrees. Requires nlat >= 2 (resp. nlon >= 2); otherwise divides by zero.
    double dlat() const { return (lat_north - lat_south) / (nlat - 1); }
    double dlon() const { return (lon_east - lon_west) / (nlon - 1); }
    // Bilinear interpolation of the height at geodetic latitude lat_deg and longitude lon_deg, both
    // in degrees (longitude may be given in any range; it is wrapped to [0, 360)). Returns NaN outside
    // [lat_south, lat_north] x [lon_west, lon_east]. Does not throw.
    double interpolate(double lat_deg, double lon_deg) const;
    // Coarser grid whose nodes are the averages over a centred window of about block_deg degrees
    // (rounded to a whole number of grid steps, at least one step). The windows are truncated at
    // the grid edges. The output nodes coincide with every k-th input node, counted from the north-
    // west corner.
    HeightGrid block_mean(double block_deg) const;
  };

  // Reads a GRAVSOFT grid. The header holds six numbers: lat1 lat2 lon1 lon2 dlat dlon. lat1 and
  // lon1 are taken as the south and west bounds, lat2 and lon2 as the north and east bounds. The
  // step is recomputed from the bounds (the file's dlat/dlon may be truncated), with
  // nlat = round((lat2 - lat1) / dlat) + 1 and similarly for longitude. Values follow row by row
  // from north to south, west to east. Throws std::runtime_error on I/O errors, an invalid header,
  // a value count that does not match the dimensions, or non-numeric values.
  HeightGrid load_gravsoft_grid(const std::string &path);

  // Abstract source of topographic heights for the computation at a set of points.
  class HeightSource
  {
  public:
    virtual ~HeightSource() = default;
    // Fills H[i] for every point i with the height at the point. Inputs:
    //   pts: evaluation points; only pts.lon (radians) is used by the grid implementation.
    //   lat_geodetic_deg: geodetic latitude of each point in DEGREES (same size as pts).
    //   o: synthesis options (ignored by the grid implementation).
    // H is resized to pts.size(). Points outside the source domain receive NaN.
    virtual void evaluate(const Points &pts, const AlignedVec<double> &lat_geodetic_deg, std::vector<double> &H, const SynthOptions &o) const = 0;
  };
  // Creates a HeightSource that interpolates `grid`. If block_deg > 0 the grid is first reduced by
  // HeightGrid::block_mean(block_deg); block_deg <= 0 uses the grid as given. Takes ownership of grid.
  std::shared_ptr<HeightSource> height_from_grid(HeightGrid grid, double block_deg);
}
