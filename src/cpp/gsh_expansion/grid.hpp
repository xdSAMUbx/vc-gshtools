// grid.hpp: synthesis on latitude-longitude grids.
// When the grid covers the full circle with uniform longitude spacing, the longitude sum is done with a 1-D
// complex-to-real FFT (FFTW). Regional grids use a direct sum with a rotation recurrence for cos/sin(mλ).
// Opposite latitudes ±φ with the same radius are computed from one set of lumps (parity of the terms in n - m).
// Units: latitudes and longitudes in radians, radii in metres, outputs dimensionless raw sums (see shs.hpp).
#pragma once
#include <string>
#include <vector>
#include "shs.hpp"

namespace gsh
{
  // A latitude-longitude grid. Rows are latitudes; columns are longitudes at lon0 + j * dlon.
  struct Grid
  {
    AlignedVec<double> lat; // spherical latitude of each row, radians; rows ordered north to south
    AlignedVec<double> r;   // one radius per row, metres
    double lon0 = 0.0;      // first longitude, radians
    double dlon = 0.0;      // longitude step, radians
    std::int64_t nlon = 0;  // number of longitudes per row
    // Number of rows (latitudes).
    std::int64_t nlat() const { return static_cast<std::int64_t>(lat.size()); }

    // Gauss-Legendre grid: L = nmax + 1 latitudes (roots of P_L, found by Newton iteration and symmetrised exactly),
    // 2L longitudes covering the full circle. Radius r for every row.
    // Same grid as CHarm (Bucha 2026).
    static Grid gauss_legendre(int nmax, double r);
    // Driscoll-Healy grid, first variant: 2L latitudes spaced by pi/(2L), starting at the north pole (lat = pi/2) with
    // the equator row set exactly to 0; 2L longitudes. Full circle. Driscoll and Healy (1994), Adv. Appl. Math. 15, 202-250.
    static Grid driscoll_healy1(int nmax, double r);
    // Driscoll-Healy grid, second variant: the same latitudes as driscoll_healy1 with 4L longitudes.
    static Grid driscoll_healy2(int nmax, double r);
    // Regional regular grid: nlat rows from lat_north to lat_south (degrees, inclusive), nlon columns from lon_west
    // to lon_east (degrees, inclusive), constant radius r (metres). The grid is a full circle only if nlon * dlon = 2π
    // (checked by full_circle()); otherwise the direct summation path is used.
    static Grid regular(double lat_north_deg, double lat_south_deg, std::int64_t nlat,
                        double lon_west_deg, double lon_east_deg, std::int64_t nlon, double r);
    // True when nlon >= 2 and nlon * dlon equals 2π within 1e-12 rad; this selects the FFT path.
    bool full_circle() const;
  };

  // Timing report of synth_grid(). lumps and fft are accumulated over threads (thread-seconds, OpenMP reduction);
  // total is wall-clock seconds.
  struct GridTiming
  {
    double lumps = 0, fft = 0, total = 0; // seconds: thread-seconds for lumps and fft, wall time for total
  };

  // Raw sums on a grid. Output: one array per field (kernels::Field order), each of nlat * nlon entries stored row by
  // row (latitude slow, longitude fast). grad = 0, 1 or 2. Rows with cos φ < 5e-3 that need derivatives are computed
  // with the full-triangle oracle (see synth_grid in grid.cpp).
  // Throws std::invalid_argument for a bad grad or a grid whose r and lat sizes differ.
  // Parallel over tiles of `lanes` rows with OpenMP; each thread has its own lump buffer and FFT spectrum.
  RawFields synth_grid(const Prepared &prep, const Grid &grid, int grad, const SynthOptions &opts = {}, GridTiming *timing = nullptr);

  // Expands a grid into scattered points, row by row (latitude slow, longitude fast). Used to validate the grid
  // path against the point path; memory and time grow as nlat * nlon, so it is slow for high degrees.
  Points grid_to_points(const Grid &grid);
}
