// ex_functionals_table.cpp: geopotential functionals at a table of points (functionals.hpp), read and written as CSV.
//
// Shows: reading a delimited table with read_table() (table_io.hpp), building Points from geodetic coordinates
// (Points::from_geodetic), requesting several functionals in one compute_functionals() call, and writing the results.
//
// Usage: ex_functionals_table <model.gfc> <points.csv> <out.csv> [nmax]
//   points.csv: one header line, then rows "id,lat_deg,lon_deg,h_m" with numeric ids and geodetic coordinates on GRS80.
//               The demo table is examples/data/points_demo.csv.
//   nmax:       highest degree, >= 2 (default: the degree of the model). Use a low value, for example 200, for a quick run.
//   The file is EGM2008 (degree 2190): loading it takes a while even for five points.
//
// Output: the CSV file, and a short table on the screen. Units are SI:
//   T        disturbing potential [m^2 s^-2]
//   zeta_ell height anomaly on the ellipsoid, T/gamma [m]
//   Dg_sa    gravity anomaly, spherical approximation [m s^-2]
//   dg       gravity disturbance, exact (|g| - gamma) [m s^-2]; undefined at the exact poles
//   V_xx, V_zz, V_xz gravitational tensor in the local frame (x north, y west, z up) [s^-2]
// N and zeta are not requested: they need a topographic height source, which this example does not use.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "functionals.hpp"
#include "table_io.hpp"

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_functionals_table: %s\n"
                         "usage: ex_functionals_table <model.gfc> <points.csv> <out.csv> [nmax]\n", why);
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
} // namespace

int main(int argc, char **argv)
{
  try
  {
    if (argc < 4 || argc > 5)
      usage("expected 3 or 4 arguments");
    const int nmax = argc > 4 ? parse_int(argv[4]) : -1;
    if (argc > 4 && nmax < 2)
      usage("nmax must be at least 2");

    const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();
    const gsh::Model model = gsh::load_gfc(argv[1], nmax);

    // Columns: 0 = id, 1 = lat_deg, 2 = lon_deg, 3 = h_m. A line with a non-numeric field is not repaired: stop instead.
    const gsh::Table tab = gsh::read_table(argv[2], ',', 1, {0, 1, 2, 3});
    if (!tab.skipped_lines.empty())
      throw std::runtime_error("non-numeric line " + std::to_string(tab.skipped_lines.front()) + " in " + argv[2]);
    const std::size_t count = tab.columns[0].size();
    if (count == 0)
      throw std::runtime_error(std::string("no points in ") + argv[2]);
    for (double lat : tab.columns[1])
      if (std::fabs(lat) >= 90.0)
        throw std::invalid_argument("dg is undefined at the exact poles; remove the pole rows from the table");

    const gsh::Points pts = gsh::Points::from_geodetic(tab.columns[1], tab.columns[2], tab.columns[3], ell.a, ell.e2);

    gsh::FunctionalRequest req;
    req.outputs = {gsh::Functional::T,    gsh::Functional::zeta_ell, gsh::Functional::Dg_sa, gsh::Functional::dg,
                   gsh::Functional::V_xx, gsh::Functional::V_zz,     gsh::Functional::V_xz};
    req.nmax = nmax;
    req.ellipsoid = ell;
    const gsh::FunctionalResult res = gsh::compute_functionals(model, pts, req);

    const auto &T = res.values.at(gsh::Functional::T);
    const auto &zeta = res.values.at(gsh::Functional::zeta_ell);
    const auto &Dg = res.values.at(gsh::Functional::Dg_sa);
    const auto &dg = res.values.at(gsh::Functional::dg);
    const auto &Vxx = res.values.at(gsh::Functional::V_xx);
    const auto &Vzz = res.values.at(gsh::Functional::V_zz);
    const auto &Vxz = res.values.at(gsh::Functional::V_xz);

    std::ofstream out(argv[3]);
    if (!out)
      throw std::runtime_error(std::string("cannot write ") + argv[3]);
    out.precision(12);
    out << "id,lat_deg,lon_deg,h_m,T_m2s2,zeta_ell_m,Dg_sa_ms2,dg_ms2,V_xx_s2,V_zz_s2,V_xz_s2\n";
    std::printf("%4s %12s %12s %10s %12s %10s %12s\n", "id", "T [m2/s2]", "zeta [m]", "Dg [mGal]", "dg [mGal]", "Vzz [E]",
                "Vxz [E]");
    for (std::size_t i = 0; i < count; ++i)
    {
      out << tab.columns[0][i] << ',' << tab.columns[1][i] << ',' << tab.columns[2][i] << ',' << tab.columns[3][i] << ','
          << T[i] << ',' << zeta[i] << ',' << Dg[i] << ',' << dg[i] << ',' << Vxx[i] << ',' << Vzz[i] << ',' << Vxz[i] << '\n';
      // Screen units: 1 mGal = 1e-5 m/s^2, 1 Eotvos = 1e-9 s^-2.
      std::printf("%4.0f %12.4f %10.4f %10.3f %12.3f %10.3f %12.3f\n", tab.columns[0][i], T[i], zeta[i], Dg[i] * 1e5,
                  dg[i] * 1e5, Vzz[i] * 1e9, Vxz[i] * 1e9);
    }
    std::printf("points=%zu, nmax=%s, polar points handled by the oracle=%lld, written to %s\n", count,
                nmax < 0 ? "model degree" : argv[4], static_cast<long long>(res.polar_points), argv[3]);
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "ex_functionals_table: %s\n", e.what());
    return 1;
  }
}
