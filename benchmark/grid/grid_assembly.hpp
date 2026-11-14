// grid_assembly.hpp: assembly of the functionals V, T and dg from the raw sums of a grid synthesis (gsh::synth_grid),
// shared by gsh_grid_bench and gsh_grid_check. The formulas are those of src/cpp/gsh_expansion/functionals.cpp
// (group V: V = GM/r S_0; group T: T = GM/r S_0 of the disturbing model, and dg from the disturbing-potential
// derivatives). Only the arithmetic is copied; no functional API is used on the grid path.
//
// Grid definitions (spherical latitudes, constant radius per row): benchmark/grid/README.md.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include "ellipsoid.hpp"
#include "functionals.hpp"
#include "grid.hpp"

namespace gridbench
{
  // Indices of the raw sums in RawFields for grad >= 1 (the order of kernels::Field, see src/cpp/gsh_legendre/kernels/kernels.hpp).
  inline constexpr int kS0 = 0; // sum q_n A_n                    (potential)
  inline constexpr int kSR = 1; // sum q_n (n+1) A_n             (radial derivative)
  inline constexpr int kSP = 2; // sum q_n A^phi_n               (latitude derivative)
  inline constexpr int kSL = 3; // sum q_n A^lambda_n            (longitude derivative, includes m)

  // V = GM/r S_0 (functionals.cpp, group V), or T = GM/r S_0 of the disturbing model (group T). Units m² s⁻².
  inline double potential(double GM, double r, double S0)
  {
    return GM / r * S0;
  }

  // Exact gravity disturbance dg = |g| - gamma (functionals.cpp, group T, field dg), in m s⁻². The disturbing gravity
  // vector delta has components (GM/r² S_P, -GM/(r² cos φ) S_L, -GM/r² S_R); gamma is the normal gravity vector of the
  // ellipsoid at the spherical latitude lat and radius r. Evaluated as (2 γ·δ + |δ|²)/(|γ+δ| + γ).
  inline double gravity_disturbance(const gsh::Ellipsoid &ell, double GM, double lat, double r, double SR, double SP, double SL)
  {
    const gsh::Ellipsoid::Normal nrm = ell.normal_gravity(lat, r);
    const double k = GM / (r * r);
    const double c = std::cos(lat);
    const double dx = k * SP;
    const double dy = -k * SL / c;
    const double dz = -k * SR;
    const double gx = nrm.g_x + dx, gz = nrm.g_z + dz;
    const double num = 2.0 * (nrm.g_x * dx + nrm.g_z * dz) + dx * dx + dy * dy + dz * dz;
    return num / (std::sqrt(gx * gx + dy * dy + gz * gz) + nrm.gamma);
  }

  // Grid of a case name: gl (Gauss-Legendre, nmax + 1 rows), dh1 and dh2 (Driscoll-Healy, 2(nmax+1) rows) or regional
  // (201 x 201, lat 45 to 35 deg north to south, lon -110 to -100 deg, step 0.05 deg). Constant radius r.
  inline gsh::Grid make_grid(const std::string &name, int nmax, double r)
  {
    if (name == "gl")
      return gsh::Grid::gauss_legendre(nmax, r);
    if (name == "dh1")
      return gsh::Grid::driscoll_healy1(nmax, r);
    if (name == "dh2")
      return gsh::Grid::driscoll_healy2(nmax, r);
    if (name == "regional")
      return gsh::Grid::regular(45.0, 35.0, 201, -110.0, -100.0, 201, r);
    throw std::invalid_argument("unknown case '" + name + "' (use gl, dh1, dh2 or regional)");
  }

  // Row indices selected for the benchmark, sorted and unique.
  //   global grids: the rows closest to 0, ±30, ±60, ±85 and ±89.5 deg, plus rows 0 and nlat - 1;
  //   regional grid: rows 0, (nlat-1)/4, (nlat-1)/2, 3(nlat-1)/4 and nlat-1 (45, 42.5, 40, 37.5 and 35 deg).
  inline std::vector<std::int64_t> auto_rows(const gsh::Grid &g, bool regional)
  {
    const std::int64_t n = g.nlat();
    std::vector<std::int64_t> rows;
    if (regional)
    {
      for (std::int64_t k = 0; k <= 4; ++k)
        rows.push_back(k * (n - 1) / 4);
    }
    else
    {
      const double targets[] = {0.0, 30.0, -30.0, 60.0, -60.0, 85.0, -85.0, 89.5, -89.5};
      for (double t : targets)
      {
        std::int64_t best = 0;
        double dbest = 1e300;
        for (std::int64_t i = 0; i < n; ++i)
        {
          const double d = std::abs(gsh::rad2deg(g.lat[i]) - t);
          if (d < dbest)
          {
            dbest = d;
            best = i;
          }
        }
        rows.push_back(best);
      }
      rows.push_back(0);
      rows.push_back(n - 1);
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    return rows;
  }

  // Parses --rows: "auto" or a comma-separated list of row indices in [0, nlat). Throws std::invalid_argument.
  inline std::vector<std::int64_t> parse_rows(const std::string &spec, const gsh::Grid &g, bool regional)
  {
    if (spec == "auto")
      return auto_rows(g, regional);
    std::vector<std::int64_t> rows;
    std::size_t pos = 0;
    while (pos <= spec.size())
    {
      const std::size_t comma = spec.find(',', pos);
      const std::string tok = spec.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
      if (tok.empty())
        throw std::invalid_argument("--rows: empty entry in '" + spec + "'");
      std::size_t used = 0;
      const long long v = std::stoll(tok, &used);
      if (used != tok.size() || v < 0 || v >= g.nlat())
        throw std::invalid_argument("--rows: '" + tok + "' is not a row index in [0, " + std::to_string(g.nlat()) + ")");
      rows.push_back(static_cast<std::int64_t>(v));
      if (comma == std::string::npos)
        break;
      pos = comma + 1;
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    return rows;
  }

  // Grid made of the given rows only, with the same longitudes and radii. Used to build Points for the selected rows.
  inline gsh::Grid sub_grid(const gsh::Grid &g, const std::vector<std::int64_t> &rows)
  {
    gsh::Grid s;
    s.lon0 = g.lon0;
    s.dlon = g.dlon;
    s.nlon = g.nlon;
    for (std::int64_t i : rows)
    {
      s.lat.push_back(g.lat[i]);
      s.r.push_back(g.r[i]);
    }
    return s;
  }
}
