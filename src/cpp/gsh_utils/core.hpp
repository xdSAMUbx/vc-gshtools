// core.hpp: shared numeric constants, angle conversions, triangular packing of spherical-harmonic
// coefficients, 64-byte aligned containers, and the geodetic-to-spherical helper.
//
// Conventions used throughout gsh:
//   - Angles given as "_deg" are degrees; internal trigonometry uses radians.
//   - Lengths are metres; GM is in m^3/s^2.
//   - Spherical-harmonic coefficients are fully normalized (unless stated otherwise) and stored
//     in a triangular layout indexed by idx_tri(n, m, nmax) (see below).
//
// Scaling of large Legendre values: the MSFC algorithm stores values as x * R^(-k) with
// R = 2^500 so that intermediate results never overflow or underflow (modified standard forward column, MSFC;
// the exponent-extension idea is that of Fukushima 2012).
#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <new>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace gsh::detail
{
  // Size in bytes of an open binary stream, 64-bit safe on MSVC (_fseeki64) and POSIX (fseeko/ftello).
  // Leaves the stream at its start.
  inline long long stream_size(std::FILE *f)
  {
#if defined(_MSC_VER)
    _fseeki64(f, 0, SEEK_END);
    const long long n = _ftelli64(f);
    _fseeki64(f, 0, SEEK_SET);
#else
    fseeko(f, 0, SEEK_END);
    const long long n = static_cast<long long>(ftello(f));
    fseeko(f, 0, SEEK_SET);
#endif
    return n;
  }
} // namespace gsh::detail

namespace gsh
{
  inline constexpr double PI = 3.14159265358979323846;

  // Degrees to radians (multiplicative, no range reduction).
  inline double deg2rad(double v) { return v * (PI / 180.0); }
  // Radians to degrees (multiplicative, no range reduction).
  inline double rad2deg(double v) { return v * (180.0 / PI); }

  // MSFC exponent rescaling: value = x * R^(-k),
  // with R = 2^RBITS. R_BIG and R_INV are the exact powers of two R and 1/R.
  inline constexpr int RBITS = 500;
  inline constexpr double R_BIG = 0x1p500;
  inline constexpr double R_INV = 0x1p-500;

  // Triangular packing of (n, m) pairs with 0 <= m <= n <= nmax, stored column by column:
  // column m holds n = m, m+1, ..., nmax contiguously. Total length for degree nmax.
  // All index helpers use 64-bit arithmetic and perform no bounds checking.
  inline constexpr std::int64_t tri_size(int nmax)
  {
    const std::int64_t N = nmax;
    return (N + 1) * (N + 2) / 2;
  }
  // Offset of the first element of column m, i.e. the element (n = m, m).
  inline constexpr std::int64_t col_offset(int m, int nmax)
  {
    const std::int64_t M = m, N = nmax;
    return M * (N + 1) - M * (M - 1) / 2;
  }
  // Index of (n, m) in the packed triangle of degree nmax. Does not check 0 <= m <= n <= nmax.
  inline constexpr std::int64_t idx_tri(int n, int m, int nmax) { return col_offset(m, nmax) + (n - m); }

  // Allocator that returns 64-byte aligned storage (one cache line, one AVX-512 ZMM register).
  // Uses _aligned_malloc on MSVC/clang-cl and std::aligned_alloc elsewhere.
  // Throws std::bad_alloc if the allocation fails.
  template <class T>
  struct AlignedAllocator
  {
    using value_type = T;
    AlignedAllocator() noexcept = default;
    template <class U>
    AlignedAllocator(const AlignedAllocator<U> &) noexcept {}
    T *allocate(std::size_t n)
    {
#if defined(_MSC_VER)
      void *p = _aligned_malloc(n * sizeof(T), 64);
#else
      void *p = std::aligned_alloc(64, ((n * sizeof(T) + 63) / 64) * 64);
#endif
      if (!p)
        throw std::bad_alloc();
      return static_cast<T *>(p);
    }
    void deallocate(T *p, std::size_t) noexcept
    {
#if defined(_MSC_VER)
      _aligned_free(p);
#else
      std::free(p);
#endif
    }
    template <class U>
    bool operator==(const AlignedAllocator<U> &) const noexcept { return true; }
  };
  // std::vector whose data pointer is 64-byte aligned.
  template <class T>
  using AlignedVec = std::vector<T, AlignedAllocator<T>>;

  // Geocentric spherical coordinates. Units: lat (geocentric latitude) and lon in RADIANS,
  // r (geocentric radius) in metres.
  struct Spherical
  {
    double lat, lon, r;
  };

  // Geodetic coordinates on the reference ellipsoid (a, e2) to geocentric spherical coordinates.
  // Inputs: lat_deg, lon_deg in degrees; h = ellipsoidal height in metres; a = semi-major axis
  // [m]; e2 = first eccentricity squared [-]. Output: Spherical with lat and lon in radians and
  // r in metres. Closed form (no iteration). Does not throw.
  // Prefer gsh::coords (coordinates.hpp) for new code; this is the shared implementation that
  // Points::from_geodetic uses.
  Spherical geodetic_to_spherical(double lat_deg, double lon_deg, double h, double a, double e2);
}
