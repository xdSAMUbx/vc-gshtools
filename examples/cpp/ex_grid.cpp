// ex_grid.cpp: synthesis on latitude-longitude grids (grid.hpp) of the potential V and the disturbing potential T.
//
// Shows: a global Gauss-Legendre grid (Grid::gauss_legendre), a regional regular grid (Grid::regular), the synthesis
// with synth_grid() and its GridTiming report, the disturbing model (disturbing_model), and the conversion of the raw
// sums S_0 into V = GM/r S_0 and T = GM/r S_0 of the disturbing model (the formula of functionals.cpp). The library has
// no public helper for this conversion on grids, so it is written out here.
//
// Usage: ex_grid <model.gfc> <regional.csv> [nmax=90]
//   nmax: truncation degree of the model and degree of the global grid, >= 2. The Gauss-Legendre grid has nmax+1
//         latitudes and 2(nmax+1) longitudes. The regional grid is fixed: 201 x 201 nodes, latitudes 45 to 35 deg N
//         and longitudes 110 to 100 deg W (step 0.05 deg). Its values are written to regional.csv, one node per line.
//
// Output: grid sizes, synthesis timing, sample values of V and T at the radius R of the model (m^2 s^-2), and the
// regional CSV with columns lat_sph_deg (spherical latitude), lon_deg, V_m2s2, T_m2s2. Timing lumps and fft are
// thread-seconds; total is wall-clock time.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include "ellipsoid.hpp"
#include "grid.hpp"
#include "model.hpp"
#include "shs.hpp"

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_grid: %s\nusage: ex_grid <model.gfc> <regional.csv> [nmax=90]\n", why);
    std::exit(2);
  }

  int parse_int(const char *s)
  {
    char *end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0')
      usage("nmax must be an integer");
    return static_cast<int>(v);
  }

  struct GridFields
  {
    gsh::AlignedVec<double> V, T; // m^2 s^-2, one entry per node, row by row (latitude slow, longitude fast)
    gsh::GridTiming tv, tt;       // timing of the synthesis of the model and of the disturbing model
  };

  // Synthesises S_0 of the model and of the disturbing model on the grid, then V = GM/r S_0 and T = GM/r S_0.
  // Each grid row has its own radius r (grid.r), so the factor GM/r is taken per row.
  GridFields synthesise(const gsh::Prepared &model_prep, const gsh::Prepared &dist_prep, const gsh::Grid &grid, double GM)
  {
    GridFields f;
    const gsh::RawFields sv = gsh::synth_grid(model_prep, grid, 0, {}, &f.tv);
    const gsh::RawFields st = gsh::synth_grid(dist_prep, grid, 0, {}, &f.tt);
    const std::int64_t nlon = grid.nlon;
    const std::size_t count = static_cast<std::size_t>(grid.nlat() * nlon);
    f.V.resize(count);
    f.T.resize(count);
    for (std::int64_t i = 0; i < grid.nlat(); ++i)
    {
      const double k = GM / grid.r[static_cast<std::size_t>(i)];
      for (std::int64_t j = 0; j < nlon; ++j)
      {
        const std::size_t idx = static_cast<std::size_t>(i * nlon + j);
        f.V[idx] = k * sv.f[0][idx];
        f.T[idx] = k * st.f[0][idx];
      }
    }
    return f;
  }

  void print_timing(const char *name, const GridFields &f)
  {
    std::printf("%s synthesis (model + disturbing model): wall %.3f s, lumps %.3f thread-s, fft %.3f thread-s\n",
                name, f.tv.total + f.tt.total, f.tv.lumps + f.tt.lumps, f.tv.fft + f.tt.fft);
  }

  void print_sample(const gsh::Grid &grid, const GridFields &f, std::int64_t row, std::int64_t col)
  {
    const std::size_t idx = static_cast<std::size_t>(row * grid.nlon + col);
    std::printf("  row %4lld (lat_sph %9.4f deg), col %4lld: V = %.6e m^2/s^2, T = %.6e m^2/s^2\n",
                static_cast<long long>(row), gsh::rad2deg(grid.lat[static_cast<std::size_t>(row)]),
                static_cast<long long>(col), f.V[idx], f.T[idx]);
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    if (argc < 3 || argc > 4)
      usage("expected 2 or 3 arguments");
    const int nmax = argc > 3 ? parse_int(argv[3]) : 90;
    if (nmax < 2)
      usage("nmax must be at least 2");

    const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();
    const gsh::Model model = gsh::load_gfc(argv[1], nmax);
    const gsh::Model dist = gsh::disturbing_model(model, ell);
    const gsh::Prepared prep_model = gsh::prepare(model);
    const gsh::Prepared prep_dist = gsh::prepare(dist);
    const double R = model.R;

    // Global grid on the sphere of radius R.
    const gsh::Grid gl = gsh::Grid::gauss_legendre(nmax, R);
    const GridFields g = synthesise(prep_model, prep_dist, gl, model.GM);
    std::printf("Gauss-Legendre grid, nmax=%d: %lld latitudes x %lld longitudes = %lld nodes\n", nmax,
                static_cast<long long>(gl.nlat()), static_cast<long long>(gl.nlon), static_cast<long long>(gl.nlat() * gl.nlon));
    print_timing("global", g);
    print_sample(gl, g, 0, 0);
    print_sample(gl, g, gl.nlat() / 2, 0);

    // Regional regular grid: direct summation path (it is not a full circle).
    const gsh::Grid reg = gsh::Grid::regular(45.0, 35.0, 201, -110.0, -100.0, 201, R);
    const GridFields q = synthesise(prep_model, prep_dist, reg, model.GM);
    std::printf("Regional grid: %lld x %lld nodes, full circle = %s\n", static_cast<long long>(reg.nlat()),
                static_cast<long long>(reg.nlon), reg.full_circle() ? "yes" : "no");
    print_timing("regional", q);
    print_sample(reg, q, 0, 0);
    print_sample(reg, q, reg.nlat() - 1, reg.nlon - 1);

    std::ofstream out(argv[2]);
    if (!out)
      throw std::runtime_error(std::string("cannot write ") + argv[2]);
    out.precision(12);
    out << "lat_sph_deg,lon_deg,V_m2s2,T_m2s2\n";
    for (std::int64_t i = 0; i < reg.nlat(); ++i)
      for (std::int64_t j = 0; j < reg.nlon; ++j)
      {
        const std::size_t idx = static_cast<std::size_t>(i * reg.nlon + j);
        out << gsh::rad2deg(reg.lat[static_cast<std::size_t>(i)]) << ','
            << gsh::rad2deg(reg.lon0 + static_cast<double>(j) * reg.dlon) << ',' << q.V[idx] << ',' << q.T[idx] << '\n';
      }
    std::printf("regional values written to %s\n", argv[2]);
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "ex_grid: %s\n", e.what());
    return 1;
  }
}
