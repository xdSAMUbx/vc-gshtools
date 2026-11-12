// gsh_meridian.cpp: gsh side of the meridional benchmark against GrafLab (see benchmark/meridian/README.md).
//
// Evaluates, at 181 points of one meridian (geodetic latitude -90..+90 deg in steps of 1 deg, one longitude,
// ellipsoidal height 0 m, GRS80), the functionals of the project API (functionals.md, compute_functionals):
//   T           disturbing potential                         [m^2 s^-2]   Functional::T
//   zeta_ell    T / gamma at the point (height anomaly ell)   [m]          Functional::zeta_ell
//   Dg_sa       gravity anomaly, spherical approximation      [m s^-2]     Functional::Dg_sa
//   dg          gravity disturbance, exact (|g| - gamma)      [m s^-2]     Functional::dg
//   V_xx,V_yy,V_zz,V_xz  gravitational tensor, LNOF frame     [s^-2]       Functional::V_xx, ...
// Units are SI in the CSV. The GrafLab side is in its own units (mGal, Eotvos); see README.md for the conversions.
//
// Two compute_functionals calls: the tensor and T/zeta_ell/Dg_sa on all 181 points, and dg on the non-polar points
// only. The API rejects dg at an exact pole (1/cos(phi) in g_y), so dg is written as "nan" at lat = +-90 deg.
//
// Usage:
//   gsh_meridian --gfc <file.gfc> --out <csv> [--nmax N] [--lon DEG]
//     --gfc   ICGEM .gfc (tide-free EGM2008 in this benchmark)
//     --nmax  maximum degree (default 2190; must not exceed the file's degree)
//     --lon   longitude of the meridian in degrees (default 0)
//     --out   result CSV, written with 17 significant digits (round-trip exact for double)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "functionals.hpp"

using namespace gsh;

namespace
{
  constexpr int kPoints = 181; // latitude -90, -89, ..., +90

  struct Args
  {
    std::string gfc, out;
    int nmax = 2190;
    double lon_deg = 0.0;
  };

  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "gsh_meridian: %s\n", why);
    std::fprintf(stderr, "usage: gsh_meridian --gfc <file.gfc> --out <csv> [--nmax N] [--lon DEG]\n");
    std::exit(2);
  }

  Args parse_args(int argc, char **argv)
  {
    Args a;
    for (int i = 1; i < argc; ++i)
    {
      const std::string key = argv[i];
      auto value = [&]() -> std::string
      {
        if (i + 1 >= argc)
          usage(("missing value after " + key).c_str());
        return argv[++i];
      };
      if (key == "--gfc")
        a.gfc = value();
      else if (key == "--out")
        a.out = value();
      else if (key == "--nmax")
        a.nmax = std::atoi(value().c_str());
      else if (key == "--lon")
        a.lon_deg = std::strtod(value().c_str(), nullptr);
      else
        usage(("unknown option " + key).c_str());
    }
    if (a.gfc.empty() || a.out.empty())
      usage("--gfc and --out are required");
    if (a.nmax < 2)
      usage("--nmax must be at least 2");
    return a;
  }

  // Result of the benchmark, one entry per point (point order = latitude -90 .. +90).
  struct Row
  {
    double lat_deg, lon_deg, h_m;
    double T, zeta_ell, Dg_sa, dg;   // SI
    double Vxx, Vyy, Vzz, Vxz;       // s^-2
  };

  void write_csv(const std::string &path, const std::vector<Row> &rows)
  {
    std::ofstream f(path);
    if (!f)
      throw std::runtime_error("cannot write " + path);
    f << "lat_deg,lon_deg,h_m,T_m2s2,zeta_ell_m,Dg_sa_ms2,dg_ms2,V_xx_s2,V_yy_s2,V_zz_s2,V_xz_s2\n";
    char buf[1024];
    for (const Row &r : rows)
    {
      // dg at the poles is not defined by the API: written as the token "nan".
      char dg_txt[64];
      if (std::isnan(r.dg))
        std::snprintf(dg_txt, sizeof dg_txt, "nan");
      else
        std::snprintf(dg_txt, sizeof dg_txt, "%.17g", r.dg);
      std::snprintf(buf, sizeof buf, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%.17g,%.17g,%.17g,%.17g\n", r.lat_deg,
                    r.lon_deg, r.h_m, r.T, r.zeta_ell, r.Dg_sa, dg_txt, r.Vxx, r.Vyy, r.Vzz, r.Vxz);
      f << buf;
    }
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    const Args a = parse_args(argc, argv);
    const Ellipsoid ell = Ellipsoid::GRS80();

    const auto t_load0 = std::chrono::steady_clock::now();
    const Model model = load_gfc(a.gfc, a.nmax);
    const double load_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_load0).count();

    std::vector<double> lat(kPoints), lon(kPoints, a.lon_deg), h(kPoints, 0.0);
    for (int i = 0; i < kPoints; ++i)
      lat[i] = -90.0 + static_cast<double>(i);
    const Points pts = Points::from_geodetic(lat, lon, h, ell.a, ell.e2);

    // ---- Request 1: all points ------------------------------------------------------------------------------
    FunctionalRequest rq;
    rq.outputs = {Functional::T, Functional::zeta_ell, Functional::Dg_sa, Functional::V_xx,
                  Functional::V_yy, Functional::V_zz, Functional::V_xz};
    rq.nmin = 0;
    rq.nmax = a.nmax;
    rq.ellipsoid = ell;
    const auto t_c1 = std::chrono::steady_clock::now();
    const FunctionalResult r1 = compute_functionals(model, pts, rq);
    const double s_c1 = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_c1).count();

    // ---- Request 2: dg on the non-polar points (the API rejects 1/cos(phi) at exact poles) -----------------
    std::vector<int> nonpolar;
    std::vector<double> lat2, lon2, h2;
    for (int i = 0; i < kPoints; ++i)
    {
      if (std::fabs(lat[i]) < 90.0)
      {
        nonpolar.push_back(i);
        lat2.push_back(lat[i]);
        lon2.push_back(lon[i]);
        h2.push_back(h[i]);
      }
    }
    const Points pts2 = Points::from_geodetic(lat2, lon2, h2, ell.a, ell.e2);
    FunctionalRequest rq2;
    rq2.outputs = {Functional::dg};
    rq2.nmin = 0;
    rq2.nmax = a.nmax;
    rq2.ellipsoid = ell;
    const auto t_c2 = std::chrono::steady_clock::now();
    const FunctionalResult r2 = compute_functionals(model, pts2, rq2);
    const double s_c2 = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_c2).count();

    // ---- Assemble and write ------------------------------------------------------------------------------------
    const auto &T = r1.values.at(Functional::T);
    const auto &zeta = r1.values.at(Functional::zeta_ell);
    const auto &dgsa = r1.values.at(Functional::Dg_sa);
    const auto &vxx = r1.values.at(Functional::V_xx);
    const auto &vyy = r1.values.at(Functional::V_yy);
    const auto &vzz = r1.values.at(Functional::V_zz);
    const auto &vxz = r1.values.at(Functional::V_xz);
    const auto &dg2 = r2.values.at(Functional::dg);

    std::vector<Row> rows(kPoints);
    for (int i = 0; i < kPoints; ++i)
    {
      Row &r = rows[i];
      r.lat_deg = lat[i];
      r.lon_deg = lon[i];
      r.h_m = h[i];
      r.T = T[i];
      r.zeta_ell = zeta[i];
      r.Dg_sa = dgsa[i];
      r.dg = std::nan("");
      r.Vxx = vxx[i];
      r.Vyy = vyy[i];
      r.Vzz = vzz[i];
      r.Vxz = vxz[i];
    }
    for (std::size_t k = 0; k < nonpolar.size(); ++k)
      rows[nonpolar[k]].dg = dg2[k];

    write_csv(a.out, rows);

    std::printf("model: %s  nmax=%d  GM=%.10g  R=%.10g  load_s=%.3f\n", a.gfc.c_str(), a.nmax, model.GM, model.R,
                load_s);
    std::printf("request 1 (T, zeta_ell, Dg_sa, V_xx, V_yy, V_zz, V_xz; %d points): %.3f s\n", kPoints, s_c1);
    std::printf("request 2 (dg; %zu non-polar points): %.3f s\n", nonpolar.size(), s_c2);
    std::printf("polar points handled by the oracle (request 1): %lld\n", static_cast<long long>(r1.polar_points));
    std::printf("wrote %s\n", a.out.c_str());
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "gsh_meridian error: %s\n", e.what());
    return 1;
  }
}
