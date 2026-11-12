// synth_types.hpp: input types of the spherical-harmonic synthesis: scattered evaluation points,
// instruction-set (ISA) selector, and the synthesis options.
//
// These types live in the gsh_utils base layer so that the HeightSource interface (height.hpp)
// does not depend on gsh_expansion. The synthesis algorithms themselves are in
// gsh_expansion/shs.hpp.
#pragma once
#include <cstdint>
#include <vector>
#include "core.hpp"

namespace gsh
{
  // Scattered evaluation points. Geocentric latitude and longitude are stored in radians and the
  // geocentric radius in metres. The three arrays always have the same length.
  struct Points
  {
    AlignedVec<double> lat, lon, r; // radians, radians, metres
    // Number of points.
    std::int64_t size() const { return static_cast<std::int64_t>(lat.size()); }
    // Builds Points from geocentric spherical coordinates: lat_deg and lon_deg in degrees, r in
    // metres. Throws std::invalid_argument if the three input vectors differ in length.
    static Points from_spherical(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg, const std::vector<double> &r);
    // Builds Points from geodetic latitude and longitude in degrees and ellipsoidal height in
    // metres, on the ellipsoid with semi-major axis a [m] and first eccentricity squared e2 [-].
    // The stored latitude is the geocentric latitude (see geodetic_to_spherical, core.hpp).
    // Throws std::invalid_argument if the three input vectors differ in length.
    static Points from_geodetic(const std::vector<double> &lat_deg, const std::vector<double> &lon_deg, const std::vector<double> &h, double a, double e2);
  };

  // Instruction-set selection for the CPU kernels. Auto picks the best ISA the CPU supports.
  enum class Isa
  {
    Auto,
    Scalar,
    Avx2,
    Avx512
  };
  // Tuning and behaviour switches for the synthesis. The defaults select the standard path.
  struct SynthOptions
  {
    int threads = 0;            // OpenMP threads; 0 = the OpenMP maximum
    int tile_vectors = 0;       // SIMD vectors per tile (1, 2, 4, 8 or 16); 0 = automatic
    Isa isa = Isa::Auto;        // ISA for the CPU kernels
    // Polar treatment: when true, the threshold max(100, 0.01 nmax) selects the polar path (shs.cpp).
    bool polar_optimization = true;
    // When true, a term is not accumulated where |P| < 2^-200 in all lanes (kernels.hpp).
    bool skip_forbidden = true;
    // When true, cos(m lambda) and sin(m lambda) are obtained by rotation with periodic resets
    // (kernels.hpp) instead of being evaluated directly for every m.
    bool rotate_trig = true;
    // When true, points are sorted by |latitude| before they are grouped into SIMD tiles, so that
    // the lanes of one tile share the same forbidden zone (shs.cpp).
    bool sort_latitude = true;
    int parallel_orders = -1;   // parallel decomposition: -1 automatic, 0 by tile, 1 by order m
  };
}
