// height.cpp: bilinear interpolation and block averaging of height grids, the GRAVSOFT reader and the
// grid-based HeightSource. See height.hpp for the interface and conventions.
#include "height.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace gsh
{
  // Bilinear interpolation on the cell that contains (lat_deg, lon_deg).
  //   fy = (lat_north - lat) / dlat   (row coordinate, 0 at the north edge)
  //   fx = (lon - lon_west) / dlon    (column coordinate, 0 at the west edge)
  // The cell indices are clamped to [0, nlat-2] x [0, nlon-2], so the upper edges use the last cell
  // with weight 1 at its far node. Weights wy, wx are the fractional offsets inside the cell.
  double HeightGrid::interpolate(double lat_deg, double lon_deg) const
  {
    double lon = std::fmod(lon_deg, 360.0);
    if (lon < 0)
      lon += 360.0;
    if (lat_deg < lat_south || lat_deg > lat_north || lon < lon_west || lon > lon_east)
      return std::numeric_limits<double>::quiet_NaN();
    const double fy = (lat_north - lat_deg) / dlat(), fx = (lon - lon_west) / dlon();
    const int iy = std::min(static_cast<int>(fy), nlat - 2), ix = std::min(static_cast<int>(fx), nlon - 2);
    const double wy = fy - iy, wx = fx - ix;
    const auto at = [&](int y, int x) -> double
    { return h[static_cast<std::size_t>(y) * nlon + x]; };
    return (1 - wy) * ((1 - wx) * at(iy, ix) + wx * at(iy, ix + 1)) +
           wy * ((1 - wx) * at(iy + 1, ix) + wx * at(iy + 1, ix + 1));
  }

  // Block averaging. The step in grid cells is k = max(1, round(block_deg / step)) for each axis.
  // Each output node (y, x) corresponds to input node (y*ky, x*kx) and averages the input nodes
  // within ky/2 rows and kx/2 columns of it (integer division), clipped to the grid. The output
  // extent is nlat' = (nlat-1)/ky + 1 nodes, and the output south/east bounds are recomputed from the
  // coarse step so that they stay on input nodes.
  HeightGrid HeightGrid::block_mean(double block_deg) const
  {
    const int ky = std::max(1, static_cast<int>(std::round(block_deg / dlat())));
    const int kx = std::max(1, static_cast<int>(std::round(block_deg / dlon())));
    HeightGrid g;
    g.nlat = (nlat - 1) / ky + 1;
    g.nlon = (nlon - 1) / kx + 1;
    g.lat_north = lat_north;
    g.lat_south = lat_north - (g.nlat - 1) * ky * dlat();
    g.lon_west = lon_west;
    g.lon_east = lon_west + (g.nlon - 1) * kx * dlon();
    g.h.resize(static_cast<std::size_t>(g.nlat) * g.nlon);
    // Centred average at each node, truncated at the edges: the nodes coincide with the output nodes.
#pragma omp parallel for schedule(static)
    for (int y = 0; y < g.nlat; ++y)
      for (int x = 0; x < g.nlon; ++x)
      {
        const int y0 = std::max(0, y * ky - ky / 2), y1 = std::min(nlat - 1, y * ky + ky / 2);
        const int x0 = std::max(0, x * kx - kx / 2), x1 = std::min(nlon - 1, x * kx + kx / 2);
        double sum = 0;
        for (int j = y0; j <= y1; ++j)
          for (int i = x0; i <= x1; ++i)
            sum += h[static_cast<std::size_t>(j) * nlon + i];
        g.h[static_cast<std::size_t>(y) * g.nlon + x] = static_cast<float>(sum / ((y1 - y0 + 1) * (x1 - x0 + 1)));
      }
    return g;
  }

  // GRAVSOFT reader. Steps:
  //   1. Read the whole file and parse the six header numbers (bounds and nominal steps).
  //   2. Recompute nlat, nlon from the bounds and the steps (rounding absorbs truncated steps).
  //   3. Count the values in parallel chunks (chunks are cut at blank characters, so no token is
  //      split), turn the counts into prefix offsets, and parse each chunk into its own slots.
  // The number of values found must equal nlat * nlon exactly; otherwise the load fails.
  HeightGrid load_gravsoft_grid(const std::string &path)
  {
    std::FILE *f = std::fopen(path.c_str(), "rb");
    if (!f)
      throw std::runtime_error("Cannot open: " + path);
    _fseeki64(f, 0, SEEK_END);
    const long long size = _ftelli64(f);
    _fseeki64(f, 0, SEEK_SET);
    std::string buf(static_cast<std::size_t>(size), '\0');
    const std::size_t got = std::fread(buf.data(), 1, buf.size(), f);
    std::fclose(f);
    if (got != buf.size())
      throw std::runtime_error("Incomplete read: " + path);

    const char *p = buf.data(), *const end = p + buf.size();
    // Skips any run of whitespace (space, tab, CR, LF).
    const auto skip = [&](const char *&q)
    { while (q < end && (*q == ' ' || *q == '\t' || *q == '\r' || *q == '\n')) ++q; };
    double head[6];
    for (double &v : head)
    {
      skip(p);
      const auto r = std::from_chars(p, end, v);
      if (r.ec != std::errc())
        throw std::runtime_error("Invalid GRAVSOFT header: " + path);
      p = r.ptr;
    }
    HeightGrid g;
    g.lat_south = head[0];
    g.lat_north = head[1];
    g.lon_west = head[2];
    g.lon_east = head[3];
    // The file's step is truncated: dimensions come from rounding, and the step is derived from the bounds.
    g.nlat = static_cast<int>(std::round((g.lat_north - g.lat_south) / head[4])) + 1;
    g.nlon = static_cast<int>(std::round((g.lon_east - g.lon_west) / head[5])) + 1;
    const std::size_t total = static_cast<std::size_t>(g.nlat) * g.nlon;
    g.h.resize(total);

    // Parallel chunks of about 4 MiB. Cut points are moved forward to the next blank so that no
    // number is split between two chunks.
    const std::int64_t body = end - p;
    const int chunks = static_cast<int>(std::max<std::int64_t>(1, body / (1 << 22)));
    std::vector<const char *> cuts(chunks + 1);
    cuts[0] = p;
    cuts[chunks] = end;
    for (int c = 1; c < chunks; ++c)
    {
      const char *q = p + body * c / chunks;
      while (q < end && *q != ' ' && *q != '\n')
        ++q;
      cuts[c] = q;
    }
    // Pass 1: count the tokens in each chunk. A token starts at a non-blank that follows a blank.
    std::vector<std::int64_t> counts(chunks + 1, 0);
#pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < chunks; ++c)
    {
      std::int64_t n = 0;
      bool in = false;
      for (const char *q = cuts[c]; q < cuts[c + 1]; ++q)
      {
        const bool blank = *q == ' ' || *q == '\t' || *q == '\r' || *q == '\n';
        if (!blank && !in)
          ++n;
        in = !blank;
      }
      counts[c + 1] = n;
    }
    // Prefix sum: counts[c] becomes the index of the first value of chunk c.
    for (int c = 0; c < chunks; ++c)
      counts[c + 1] += counts[c];
    if (counts[chunks] != static_cast<std::int64_t>(total))
      throw std::runtime_error("GRAVSOFT grid: unexpected number of values in " + path);
    // Pass 2: parse each chunk into its own slots. `bad` counts tokens that are not numbers.
    int bad = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : bad)
    for (int c = 0; c < chunks; ++c)
    {
      std::int64_t idx = counts[c];
      const char *q = cuts[c];
      const char *stop = cuts[c + 1];
      while (true)
      {
        while (q < stop && (*q == ' ' || *q == '\t' || *q == '\r' || *q == '\n'))
          ++q;
        if (q >= stop)
          break;
        double v = 0;
        const auto r = std::from_chars(q, stop, v);
        if (r.ec != std::errc())
        {
          ++bad;
          break;
        }
        g.h[idx++] = static_cast<float>(v);
        q = r.ptr;
      }
    }
    if (bad)
      throw std::runtime_error("GRAVSOFT grid with non-numeric values: " + path);
    return g;
  }

  namespace
  {
    // HeightSource backed by a (possibly block-averaged) grid. Each point is interpolated bilinearly;
    // the longitude is converted from radians (Points) to degrees for the grid lookup.
    class GridSource final : public HeightSource
    {
    public:
      explicit GridSource(HeightGrid g) : grid_(std::move(g)) {}
      // The options are ignored: interpolation is a plain loop, parallelised with OpenMP.
      void evaluate(const Points &pts, const AlignedVec<double> &lat_geod, std::vector<double> &H, const SynthOptions &) const override
      {
        const std::int64_t n = pts.size();
        H.resize(n);
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < n; ++i)
          H[i] = grid_.interpolate(lat_geod[i], rad2deg(pts.lon[i]));
      }

    private:
      HeightGrid grid_;
    };
  }

  // Wraps the grid in a GridSource, after block averaging if block_deg > 0.
  std::shared_ptr<HeightSource> height_from_grid(HeightGrid grid, double block_deg)
  {
    return std::make_shared<GridSource>(block_deg > 0 ? grid.block_mean(block_deg) : std::move(grid));
  }
}
