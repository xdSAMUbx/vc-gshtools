# height: topography for N and ζ

Header `src/cpp/gsh_utils/height.hpp`, implementation `src/cpp/gsh_utils/height.cpp`. Specification §7 (internal study notes, not distributed).

The specification defines H with spherical-harmonic coefficients. The Colorado DEM is a regional grid that cannot be expanded, so there are two sources with the same interface.

## `struct HeightGrid`

| Field | Meaning |
|---|---|
| `lat_south`, `lat_north`, `lon_west`, `lon_east` | Limits in degrees. Longitudes 0..360 |
| `nlat`, `nlon` | Dimensions |
| `h` | Heights in metres, `float`, rows from north to south |

`dlat()` and `dlon()` give the step, derived from the limits and the dimensions.

**Example**

```cpp
#include "height.hpp"
#include <cstdio>
using namespace gsh;

HeightGrid g;
g.lat_south = 0;  g.lat_north = 2;   // 2 degrees of latitude
g.lon_west = 0;   g.lon_east = 2;    // 2 degrees of longitude
g.nlat = 3;       g.nlon = 3;        // 3 x 3 nodes, 1 degree apart
g.h = {0, 10, 20,  100, 110, 120,  200, 210, 220};  // row 0 is lat 2 (north)
std::printf("dlat=%.1f dlon=%.1f h(0,0)=%.0f\n", g.dlat(), g.dlon(), g.h[0]);
// Output: dlat=1.0 dlon=1.0 h(0,0)=0
```

### `double interpolate(double lat_deg, double lon_deg) const`
Bilinear interpolation. Accepts longitudes in any range and returns **NaN** outside the grid.

**Example**

```cpp
#include "height.hpp"
#include <cstdio>
using namespace gsh;

// Grid with h = 100*(2 - lat) + 10*lon, so bilinear interpolation is exact.
HeightGrid g;
g.lat_south = 0; g.lat_north = 2; g.lon_west = 0; g.lon_east = 2;
g.nlat = 3; g.nlon = 3;
g.h = {0, 10, 20, 100, 110, 120, 200, 210, 220};
std::printf("%.1f\n", g.interpolate(1.5, 0.5));    // 55.0   (between nodes)
std::printf("%.1f\n", g.interpolate(0.0, 2.0));    // 220.0  (south-east corner)
std::printf("%.1f\n", g.interpolate(1.5, 360.5));  // 55.0   (longitude wrapped to 0.5)
std::printf("%.1f\n", g.interpolate(2.5, 1.0));    // nan    (north of lat_north)
```

### `HeightGrid block_mean(double block_deg) const`
Block average over `block_deg` degrees, centred on each output node and truncated at the edges.

**Example**

```cpp
#include "height.hpp"
#include <cstdio>
using namespace gsh;

HeightGrid g;   // same 3 x 3 grid as the interpolate example
g.lat_south = 0; g.lat_north = 2; g.lon_west = 0; g.lon_east = 2;
g.nlat = 3; g.nlon = 3;
g.h = {0, 10, 20, 100, 110, 120, 200, 210, 220};
HeightGrid c = g.block_mean(2.0);   // 2 x 2 nodes, 2 degrees apart
std::printf("nlat=%d nlon=%d dlat=%.1f\n", c.nlat, c.nlon, c.dlat());  // nlat=2 nlon=2 dlat=2.0
std::printf("%.1f %.1f %.1f %.1f\n", c.h[0], c.h[1], c.h[2], c.h[3]); // 55.0 65.0 155.0 165.0
std::printf("%.1f\n", c.interpolate(1.0, 1.0));                        // 110.0 (centre)
```

## `HeightGrid load_gravsoft_grid(const std::string &path)`

Reads the GRAVSOFT format: header `lat1 lat2 lon1 lon2 dlat dlon`, and values row by row from north to south. The Colorado DEM (6001 × 9601) takes 0.2 s.

- **Dimensions:** obtained by rounding with the step in the file, and the real step is recomputed from the limits, because a truncated `0.00083333` accumulates about 2 m of error.
- **Parallelism:** counts values per chunk to fix the offset of each one, then converts.
- **Errors:** `runtime_error` if the file cannot be opened, the header is invalid, the number of values does not match, or there are non-numeric values.

**Example**

```cpp
#include "height.hpp"
#include <cstdio>
using namespace gsh;

HeightGrid dem;
try
{
  dem = load_gravsoft_grid("data/cache/tablas/Colorado_dem/Colorado_dem.grd");  // needs the DEM file; not run here
  std::printf("nlat=%d nlon=%d dlat=%.10f deg\n", dem.nlat, dem.nlon, dem.dlat());
  // Expected for the Colorado DEM: nlat=6001 nlon=9601.
}
catch (const std::runtime_error &e)   // also thrown for a missing file or a bad header
{
  std::fprintf(stderr, "DEM error: %s\n", e.what());
}
```

## `class HeightSource`

Abstract interface.

### `virtual void evaluate(pts, lat_geodetic_deg, H, opts) const`
Fills `H` with one height per point, NaN where there is no data.

- `pts`: evaluation points (spherical coordinates in radians). The grid implementation uses only `pts.lon`.
- `lat_geodetic_deg`: geodetic latitude of each point in degrees, same size as `pts`. This is the latitude the grid implementation looks up.
- `opts`: synthesis options (ignored by the grid implementation).

**Example** (a minimal subclass that implements `evaluate`)

```cpp
#include "height.hpp"
#include <cmath>
#include <limits>
using namespace gsh;

// Custom source: 500 m for |lat| <= 45 deg, NaN elsewhere.
class BandSource final : public HeightSource
{
public:
  void evaluate(const Points &pts, const AlignedVec<double> &lat, std::vector<double> &H, const SynthOptions &) const override
  {
    H.resize(static_cast<std::size_t>(pts.size()));
    for (std::size_t i = 0; i < H.size(); ++i)
      H[i] = std::fabs(lat[i]) <= 45.0 ? 500.0 : std::numeric_limits<double>::quiet_NaN();
  }
};

BandSource src;
Points pts = Points::from_geodetic({10.0, 60.0}, {0.0, 0.0}, {0.0, 0.0}, 6378137.0, 0.00669437999014);
AlignedVec<double> lat_geod{10.0, 60.0};   // geodetic degrees, same size as pts
std::vector<double> H;
src.evaluate(pts, lat_geod, H, SynthOptions{});  // H = {500.0, NaN}
```

### `std::shared_ptr<HeightSource> height_from_grid(HeightGrid grid, double block_deg)`
Source from a grid (`GridSource`, internal). With `block_deg <= 0` the grid is used without averaging. The example uses 5′ (`5.0/60.0`), consistent with the resolution of EGM2008.

**Example**

```cpp
#include "height.hpp"
#include <utility>
using namespace gsh;

HeightGrid g;   // 3 x 3 grid of the interpolate example
g.lat_south = 0; g.lat_north = 2; g.lon_west = 0; g.lon_east = 2;
g.nlat = 3; g.nlon = 3;
g.h = {0, 10, 20, 100, 110, 120, 200, 210, 220};
auto src = height_from_grid(g, 0.0);                // copies g; no averaging
auto coarse = height_from_grid(std::move(g), 2.0);  // block-averaged to 2 degrees
```

```cpp
Points pts = Points::from_geodetic({1.5, 0.0, 2.5}, {0.5, 2.0, 1.0}, {0.0, 0.0, 0.0}, 6378137.0, 0.00669437999014);
AlignedVec<double> lat_geod{1.5, 0.0, 2.5};   // geodetic degrees
std::vector<double> H;
src->evaluate(pts, lat_geod, H, SynthOptions{});  // H = {55.0, 220.0, NaN}
```

### `std::shared_ptr<HeightSource> height_from_model(Model model, int nmax = -1)`
Declared in `src/cpp/gsh_expansion/height_model.hpp` (implementation `height_model.cpp`). Source from an HC/HS model (`ModelSource`, internal). It is evaluated with `q = 1` and without the GM/r factor.

**Example**

```cpp
#include "height_model.hpp"   // height_from_model (pulls in height.hpp)
#include "model.hpp"          // Model, Model::zeros, load_gfc
using namespace gsh;

// With a real HC/HS topography model (ICGEM .gfc file; path is a placeholder):
//   auto src = height_from_model(load_gfc("path/to/topo.gfc", 2190), 2190);
// Runnable check: all coefficients zero, so H = 0 m everywhere.
auto src = height_from_model(Model::zeros(3.986004418e14, 6378137.0, 360), 360);
Points pts = Points::from_geodetic({45.0}, {100.0}, {0.0}, 6378137.0, 0.00669437999014);
AlignedVec<double> lat_geod{45.0};
std::vector<double> H;
src->evaluate(pts, lat_geod, H, SynthOptions{});  // H = {0.0}; this source ignores lat_geod and uses pts.lat
```

## Unconfirmed assumptions

H of the DEM is taken as **orthometric**, because the file does not declare it. The term `2πGρH²/γ` reaches about 2.2 m, so the source and the resolution matter: changing from 5′ to 3′ changes N by 0.061 m rms, with extremes of up to 0.76 m.

## Putting it together

Typical workflow: load the regional DEM, reduce it to 5′ blocks, and evaluate heights at geodetic points. The result feeds the N and ζ functionals (see `examples/cpp/ex_aerogravimetria.cpp`, where `rq.height = height_from_grid(dem, 5.0 / 60.0)`).

```cpp
#include "height.hpp"
#include <cstdio>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace gsh;

int main()
{
  // 1. Load the regional DEM (GRAVSOFT, heights in metres, taken as orthometric).
  HeightGrid dem;
  try { dem = load_gravsoft_grid("data/cache/tablas/Colorado_dem/Colorado_dem.grd"); }
  catch (const std::runtime_error &e) { std::fprintf(stderr, "%s\n", e.what()); return 1; }

  // 2. Block-average to 5 arc-minutes (EGM2008 resolution) and wrap as a HeightSource.
  auto src = height_from_grid(std::move(dem), 5.0 / 60.0);

  // 3. Evaluate at geodetic points. Coordinates are illustrative; pick points inside the DEM.
  const std::vector<double> lat = {39.0, 39.5}, lon = {254.0, 254.5}, h = {0.0, 0.0};
  Points pts = Points::from_geodetic(lat, lon, h, 6378137.0, 0.00669437999014);
  AlignedVec<double> lat_geod(lat.begin(), lat.end());   // the grid needs geodetic degrees
  std::vector<double> H;
  src->evaluate(pts, lat_geod, H, SynthOptions{});

  // 4. H is in metres; NaN means the point is outside the DEM.
  for (std::size_t i = 0; i < H.size(); ++i)
    std::printf("point %zu: H = %.2f m\n", i, H[i]);
  return 0;
}
```
