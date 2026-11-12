// gsh_shs_bench.cpp: driver of the gsh engine for the spherical-harmonic synthesis benchmark.
//
// Definitions: benchmark/shs/results/csv.md. The driver evaluates, at every point of a CSV set, three functionals in ONE
// pass (one run = one result CSV = one timings row per repetition):
//   V_m2s2   gravitational potential of the model V = GM/r * S0(C, S)                    [m^2 s^-2]
//   T_m2s2   disturbing potential T = GM/r * S0(dC, dS) = V - U                          [m^2 s^-2]
//   dg_mgal  exact gravity disturbance dg = |g| - gamma(r, phi)                          [mGal]
// where g is the gravity vector of the model (C, S) in the LNOF frame, with the centrifugal term, and gamma is the
// normal gravity of the GRS80 field series (Ellipsoid::normal_gravity). The definitions are those of the project
// (src/cpp/gsh_expansion/functionals.cpp, functions compute_functionals, groups V and T). The dg value in m s^-2 is
// multiplied by 1e5 to give mGal at write time.
//
// Recurrence modes (--recurrence):
//   fast  the production path: gsh::compute_functionals() with the outputs {V, T, dg} in one call (SIMD kernels,
//         OpenMP over points). The model is synthesised once with gradient order 1 (its S0 gives V, its gradient
//         sums give dg) and the disturbing model once with order 0 (its S0 gives T). No functional is computed twice.
//   sfc   per point: the fully normalised Legendre triangle P and its latitude derivative dP come from
//         legendre::sfc() (+ legendre::dphi). The series are summed in this driver with the same loop and formulas
//         as the project's oracle shs::reference_point() (src/cpp/gsh_expansion/reference.cpp).
//   msfc  as sfc, with legendre::msfc().
//   era   as sfc, with legendre::era() (extended-range arithmetic; slow, precision reference).
//
// The geodetic -> geocentric spherical conversion is done ONCE for all modes, by Points::from_geodetic (GRS80), so
// that the four recurrences evaluate the same geometry. The three functionals are identical for all modes.
//
// Timing: the timed block covers everything after the model and the points are loaded: the disturbing model,
// the recurrence tables, the per-point synthesis and the functional formulas. Wall time is std::chrono::steady_clock;
// CPU time is the process user+kernel time (GetProcessTimes on Windows, CLOCK_PROCESS_CPUTIME_ID elsewhere), so it
// includes all OpenMP threads. Load time (model and points) is written to the `notes` field as load_s=...
//
// Usage:
//   gsh_shs_bench --gfc <file.gfc> --points <csv> --recurrence sfc|msfc|era|fast --nmax N --out <csv>
//                 --timings results/timings.csv --set colorado|latsweep|grid --threads K --repeat R
// threads = 0 means the OpenMP default. The result CSV (id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal) is written
// from the last repetition, with a companion <out>.md; one timings row is appended per repetition (the header is
// written when the file does not exist yet).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <ctime>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "core.hpp"
#include "ellipsoid.hpp"
#include "functionals.hpp"
#include "legendre.hpp"
#include "model.hpp"
#include "synth_types.hpp"

using namespace gsh;

namespace
{
  // ------------------------------------------------------------------------------------------------------------
  // Command line
  // ------------------------------------------------------------------------------------------------------------
  struct Args
  {
    std::string gfc, points, recurrence, out, timings, set;
    int nmax = -1;
    int threads = 0;   // 0 = OpenMP default
    int repeat = 1;
  };

  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "gsh_shs_bench: %s\n", why);
    std::fprintf(stderr,
                 "usage: gsh_shs_bench --gfc <file.gfc> --points <csv> --recurrence sfc|msfc|era|fast --nmax N\n"
                 "                     --out <csv> --timings <csv> --set colorado|latsweep|grid --threads K --repeat R\n");
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
      else if (key == "--points")
        a.points = value();
      else if (key == "--recurrence")
        a.recurrence = value();
      else if (key == "--nmax")
        a.nmax = std::atoi(value().c_str());
      else if (key == "--out")
        a.out = value();
      else if (key == "--timings")
        a.timings = value();
      else if (key == "--set")
        a.set = value();
      else if (key == "--threads")
        a.threads = std::atoi(value().c_str());
      else if (key == "--repeat")
        a.repeat = std::atoi(value().c_str());
      else
        usage(("unknown option " + key).c_str());
    }
    if (a.gfc.empty() || a.points.empty() || a.out.empty() || a.timings.empty())
      usage("--gfc, --points, --out and --timings are required");
    if (a.recurrence != "sfc" && a.recurrence != "msfc" && a.recurrence != "era" && a.recurrence != "fast")
      usage("--recurrence must be sfc, msfc, era or fast");
    if (a.set != "colorado" && a.set != "latsweep" && a.set != "grid")
      usage("--set must be colorado, latsweep or grid");
    if (a.nmax < 0)
      usage("--nmax must be a non-negative integer");
    if (a.repeat < 1 || a.threads < 0)
      usage("--repeat must be >= 1 and --threads >= 0");
    return a;
  }

  // ------------------------------------------------------------------------------------------------------------
  // Points file: header id,lat_deg,lon_deg,h_m; one row per point (see results/csv.md).
  // ------------------------------------------------------------------------------------------------------------
  struct PointFile
  {
    std::vector<long long> id;
    std::vector<double> lat_deg, lon_deg, h_m;
  };

  PointFile read_points(const std::string &path)
  {
    std::ifstream in(path);
    if (!in)
      throw std::runtime_error("cannot open points file " + path);
    std::string line;
    if (!std::getline(in, line))
      throw std::runtime_error("empty points file " + path);
    while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
      line.pop_back();
    if (line != "id,lat_deg,lon_deg,h_m")
      throw std::runtime_error("unexpected header in " + path + ": " + line);

    PointFile pf;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r')
        line.pop_back();
      if (line.empty())
        continue;
      // Four comma-separated fields. strtod is used (not iostream) so that the parse is locale-independent
      // with the C locale, which is the default of a fresh process.
      const char *p = line.c_str();
      char *end = nullptr;
      const long long id = std::strtoll(p, &end, 10);
      if (*end != ',')
        throw std::runtime_error("bad row in " + path + ": " + line);
      const double lat = std::strtod(end + 1, &end);
      if (*end != ',')
        throw std::runtime_error("bad row in " + path + ": " + line);
      const double lon = std::strtod(end + 1, &end);
      if (*end != ',')
        throw std::runtime_error("bad row in " + path + ": " + line);
      const double h = std::strtod(end + 1, &end);
      pf.id.push_back(id);
      pf.lat_deg.push_back(lat);
      pf.lon_deg.push_back(lon);
      pf.h_m.push_back(h);
    }
    return pf;
  }

  // ------------------------------------------------------------------------------------------------------------
  // Timers
  // ------------------------------------------------------------------------------------------------------------
  double wall_now()
  {
    using Clock = std::chrono::steady_clock;
    static const Clock::time_point origin = Clock::now();
    return std::chrono::duration<double>(Clock::now() - origin).count();
  }

  // Process CPU time (user + kernel), all threads included.
  double cpu_now()
  {
#ifdef _WIN32
    FILETIME creation, exit_t, kernel, user;
    GetProcessTimes(GetCurrentProcess(), &creation, &exit_t, &kernel, &user);
    const auto to_s = [](const FILETIME &f)
    {
      ULARGE_INTEGER v;
      v.LowPart = f.dwLowDateTime;
      v.HighPart = f.dwHighDateTime;
      return static_cast<double>(v.QuadPart) * 1e-7; // 100 ns units
    };
    return to_s(kernel) + to_s(user);
#else
    timespec ts{};
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
#endif
  }

  int omp_threads_now()
  {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  }

  // ------------------------------------------------------------------------------------------------------------
  // Results: V and T in m^2 s^-2, dg in m s^-2 (converted to mGal when written).
  // ------------------------------------------------------------------------------------------------------------
  struct Outputs
  {
    std::vector<double> V, T, dg;
  };

  // Raw sums of spec section 2 at one point, the same loop as shs::reference_point() (reference.cpp):
  //   S0 = sum_n q_n A_n,  SR = sum_n q_n (n+1) A_n,  SP = sum_n q_n A^phi_n,  SL = sum_n q_n A^lambda_n,
  // with A = (C cos m lam + S sin m lam) P, A^phi = (...) dP/dphi, A^lambda = m (S cos m lam - C sin m lam) P,
  // q_n = (R/r)^n. S0 gives V (model) and T (disturbing model); the derivative sums give dg.
  struct RawSums
  {
    double S0 = 0.0, SR = 0.0, SP = 0.0, SL = 0.0;
  };

  RawSums raw_sums(const Model &model, int nmax, double lon, double r, const double *P, const double *dP,
                   bool with_derivatives, std::vector<double> &qn)
  {
    // Powers of q = R/r, same construction as reference_point().
    const double q = model.R / r;
    qn[0] = 1.0;
    for (int n = 1; n <= nmax; ++n)
      qn[n] = qn[n - 1] * q;

    // Neumaier compensated accumulators, as in reference_point(): a plain double sum over 2.4e6 terms loses about
    // sqrt(N) eps (1e-5 m²/s² in V at degree 2190). Requires /fp:precise (no fast-math).
    struct Neumaier
    {
      double s = 0.0, c = 0.0;
      void add(double x)
      {
        const double t = s + x;
        c += (std::fabs(s) >= std::fabs(x)) ? (s - t) + x : (x - t) + s;
        s = t;
      }
      double value() const { return s + c; }
    };
    Neumaier k0, kr, kp, kl;
    for (int m = 0; m <= nmax; ++m)
    {
      const double cm = std::cos(m * lon), sm = std::sin(m * lon);
      const double M = m;
      for (int n = m; n <= nmax; ++n)
      {
        const double N = n;
        const double C = model.c(n, m), S = model.s(n, m);
        const auto i = idx_tri(n, m, nmax);
        const double base = C * cm + S * sm, cross = S * cm - C * sm;
        const double A = base * P[i];
        const double Q = qn[n];
        k0.add(Q * A);
        if (with_derivatives)
        {
          const double Ap = base * dP[i];
          const double Al = M * cross * P[i];
          kr.add(Q * (N + 1) * A);
          kp.add(Q * Ap);
          kl.add(Q * Al);
        }
      }
    }
    RawSums s;
    s.S0 = k0.value();
    s.SR = kr.value();
    s.SP = kp.value();
    s.SL = kl.value();
    return s;
  }

  // Legendre triangle for the per-point modes. Same pole convention as reference_point(): at a pole t = +-1, u = 0.
  using LegendreFn = void (*)(double, double, const legendre::Coefs &, std::span<double>);

  // Per-point synthesis (sfc, msfc, era), OpenMP over points. Definitions identical to compute_functionals().
  Outputs run_per_point(const Model &model, const Model &dist, const Points &pts, int nmax, const Ellipsoid &ell,
                        LegendreFn legendre_fn, const legendre::Coefs &coefs, int threads)
  {
    const int np = static_cast<int>(pts.size());
    Outputs out;
    out.V.assign(np, 0.0);
    out.T.assign(np, 0.0);
    out.dg.assign(np, 0.0);
    const std::size_t tri = static_cast<std::size_t>(tri_size(nmax));
    const double GM = model.GM;
    const double w2 = ell.omega * ell.omega;

    if (threads > 0)
    {
#ifdef _OPENMP
      omp_set_num_threads(threads);
#endif
    }

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < np; ++i)
    {
      // Per-thread work buffers: allocated once per thread, reused for all points of that thread.
      thread_local std::vector<double> P, dP, qn;
      if (P.size() != tri)
      {
        P.assign(tri, 0.0);
        dP.assign(tri, 0.0);
      }
      if (qn.size() != static_cast<std::size_t>(nmax) + 1)
        qn.assign(static_cast<std::size_t>(nmax) + 1, 0.0);

      const double lat = pts.lat[i], lon = pts.lon[i], r = pts.r[i];
      const bool pole = std::fabs(lat) == PI / 2;
      const double t = pole ? (lat > 0 ? 1.0 : -1.0) : std::sin(lat);
      const double u = pole ? 0.0 : std::cos(lat);

      legendre_fn(t, u, coefs, std::span<double>(P.data(), P.size()));
      legendre::dphi(nmax, std::span<const double>(P.data(), P.size()), std::span<double>(dP.data(), dP.size()));

      // Model: S0 gives V (functionals.cpp, group V, PDF 1); the derivative sums give dg. One loop, both results.
      const RawSums sm = raw_sums(model, nmax, lon, r, P.data(), dP.data(), true, qn);
      out.V[i] = GM / r * sm.S0;

      // T = GM/r S0 of the disturbing model (functionals.cpp, group T, PDF 49).
      const RawSums sd = raw_sums(dist, nmax, lon, r, P.data(), dP.data(), false, qn);
      out.T[i] = GM / r * sd.S0;

      // dg = |g| - gamma (functionals.cpp, group V, with the same LNOF gravity vector).
      const double cosl = std::cos(lat), sinl = std::sin(lat);
      const double dVr = -GM / (r * r) * sm.SR;
      const double dVp = GM / r * sm.SP;
      const double dVl = GM / r * sm.SL;
      const double gx = (dVp - w2 * r * r * cosl * sinl) / r;
      const double gy = -dVl / (r * cosl);
      const double gz = dVr + w2 * r * cosl * cosl;
      const double gmod = std::sqrt(gx * gx + gy * gy + gz * gz);
      out.dg[i] = gmod - ell.normal_gravity(lat, r).gamma;
    }
    return out;
  }

  // Production path: gsh::compute_functionals with V, T and dg in one request. The library plans the passes:
  // the model at gradient order 1 (its S0 is V, its gradient sums are dg) and the disturbing model at order 0 (T).
  Outputs run_fast(const Model &model, const Points &pts, int nmax, int threads)
  {
    FunctionalRequest rq;
    rq.outputs = {Functional::V, Functional::T, Functional::dg};
    rq.nmax = nmax;
    rq.nmin = 0;   // explicit: degree 0 (C00 of the file) is part of V
    rq.ellipsoid = Ellipsoid::GRS80();
    rq.synth.threads = threads;
    const FunctionalResult res = compute_functionals(model, pts, rq);
    const auto &v = res.values.at(Functional::V);
    const auto &t = res.values.at(Functional::T);
    const auto &g = res.values.at(Functional::dg);
    Outputs out;
    out.V.assign(v.begin(), v.end());
    out.T.assign(t.begin(), t.end());
    out.dg.assign(g.begin(), g.end());
    return out;
  }

  // ------------------------------------------------------------------------------------------------------------
  // Output
  // ------------------------------------------------------------------------------------------------------------
  void write_result(const std::string &path, const PointFile &pf, const Outputs &o)
  {
    std::ofstream f(path);
    if (!f)
      throw std::runtime_error("cannot write " + path);
    f << "id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal\n";
    char buf[512];
    for (std::size_t i = 0; i < pf.id.size(); ++i)
    {
      // 15 significant digits, as in results/csv.md. dg is written in mGal (1 m/s^2 = 1e5 mGal).
      std::snprintf(buf, sizeof buf, "%lld,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n", pf.id[i], pf.lat_deg[i],
                    pf.lon_deg[i], pf.h_m[i], o.V[i], o.T[i], o.dg[i] * 1e5);
      f << buf;
    }
  }

  // Companion .md of the result: column table, exact command and the formulas used.
  void write_companion_md(const std::string &csv_path, const Args &a, const std::string &command, int n_points,
                          int threads, double load_s, const Model &model, const std::string &model_name)
  {
    std::string md = csv_path;
    if (md.size() >= 4 && md.compare(md.size() - 4, 4, ".csv") == 0)
      md.resize(md.size() - 4);
    md += ".md";
    std::ofstream f(md);
    if (!f)
      throw std::runtime_error("cannot write " + md);
    const std::size_t sep = csv_path.find_last_of("/\\");
    f << "# " << (sep == std::string::npos ? csv_path : csv_path.substr(sep + 1)) << "\n\n";
    f << "Engine: **gsh** (C++ driver `gsh_shs_bench`), recurrence `" << a.recurrence << "`, point set `" << a.set
      << "`, functionals V, T and dg computed in one pass.\n\n";
    f << "## Columns\n\n| column | unit | meaning |\n|---|---|---|\n"
         "| id | - | point identifier from the input CSV |\n"
         "| lat_deg | deg | geodetic latitude on GRS80 (input, echoed) |\n"
         "| lon_deg | deg | longitude (input, echoed) |\n"
         "| h_m | m | ellipsoidal height on GRS80 (input, echoed) |\n"
         "| V_m2s2 | m^2/s^2 | gravitational potential V = GM/r sum_n (R/r)^n sum_m (C cos m lambda + S sin m lambda) P_nm(sin phi') "
         "over n = 0..nmax, with the coefficients of the model as read from the file (C00 included); "
         "geocentric r and phi' from the GRS80 geodetic input; no normal potential subtracted, so V is not T |\n"
         "| T_m2s2 | m^2/s^2 | disturbing potential T = V - U at the point (V: model, U: GRS80 normal potential, closed "
         "form; T = GM/r sum_n (R/r)^n sum_m (dC cos m lambda + dS sin m lambda) P_nm, with the disturbing "
         "coefficients dC, dS of the project) |\n"
         "| dg_mgal | mGal | gravity disturbance dg = \\|g\\| - gamma_SH at the point (g: gravity vector of the model in the "
         "LNOF frame, including the centrifugal term; gamma_SH: normal gravity of the GRS80 zonal series of the project) |\n\n";
    f << "## Command and parameters\n\n```\n" << command << "\n```\n\n";
    f << "- recurrence: `" << a.recurrence << "`\n";
    f << "- nmax: `" << a.nmax << "`, nmin: `0`\n";
    f << "- threads: `" << threads << "` (0 on the command line means the OpenMP default)\n";
    f << "- repeat: `" << a.repeat << "` (the result is written from the last repetition)\n";
    f << "- points: `" << n_points << "`\n";
    f << "- method: " << (a.recurrence == "fast"
                             ? "`gsh::compute_functionals` with the outputs `{V, T, dg}` in one call: the model is synthesised once "
                               "with gradient order 1 (its S0 gives V, its gradient sums give dg) and the disturbing model once "
                               "with order 0 (its S0 gives T). SIMD kernels, OpenMP over points."
                             : "per point: `legendre::" + a.recurrence + "` gives the normalised triangle P_nm and "
                               "`legendre::dphi` its latitude derivative; the driver forms the raw sums of the model "
                               "(S0 for V, derivative sums for dg) and of the disturbing model (S0 for T) in one loop per model.")
      << "\n\n";
    f << "## Model\n\n";
    f << "- file: `" << a.gfc << "`\n";
    f << "- modelname: `" << model_name << "`\n";
    f << "- GM [m^3/s^2]: `" << model.GM << "`, R [m]: `" << model.R << "`\n";
    f << "- tide_system (header): `" << model.tide_system << "` (no tide conversion applied by the driver)\n";
    f << "- load_s: `" << load_s << "`\n\n";
    f << "## Notes\n\n"
         "- V is the full gravitational potential of the model. It is NOT the disturbing potential T = V - U, and no "
         "normal potential is subtracted. Compare it with GrafLab Functional 11 with nmin = 0.\n"
         "- Degree 0 is included in V: the term GM/r * C00 uses C00 from the file (1.0 in the tide-free model).\n"
         "- Geodetic to geocentric conversion: `Points::from_geodetic` (GRS80), the same for all recurrences.\n";
  }

  void append_timing(const std::string &path, const Args &a, int n_points, int threads, int repeat, double wall,
                     double cpu, double load_s, const std::string &model_name)
  {
    const bool exists = [&]()
    {
      std::ifstream test(path);
      return test.good();
    }();
    std::ofstream f(path, std::ios::app);
    if (!f)
      throw std::runtime_error("cannot append to " + path);
    if (!exists)
      f << "engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes\n";
    char buf[1024];
    // notes is one field: no commas inside it.
    std::snprintf(buf, sizeof buf, "gsh,%s,%s,%d,%d,%d,%d,%.6f,%.6f,load_s=%.3f;model=%s;out=%s;functional=V_T_dg\n",
                  a.recurrence.c_str(), a.set.c_str(), a.nmax, n_points, threads, repeat, wall, cpu, load_s,
                  model_name.c_str(), a.out.c_str());
    f << buf;
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    const Args a = parse_args(argc, argv);
    const Ellipsoid ell = Ellipsoid::GRS80();
    std::string command = "gsh_shs_bench";
    for (int i = 1; i < argc; ++i)
      command += std::string(" ") + argv[i];

    // Sfc beyond 1800: the plain-double sectorial term underflows (legendre.hpp). No guard; only a warning.
    if (a.recurrence == "sfc" && a.nmax > 1800)
      std::fprintf(stderr, "WARNING: --recurrence sfc with nmax=%d > 1800 is expected to be inaccurate "
                           "(sectorial underflow, legendre.md section 5)\n", a.nmax);

    // ---- Load: model, points, geocentric coordinates (not timed) ----------------------------------------
    const double t_load0 = wall_now();
    const Model model = load_gfc(a.gfc, a.nmax);
    const PointFile pf = read_points(a.points);
    const int np = static_cast<int>(pf.id.size());
    const Points pts = Points::from_geodetic(pf.lat_deg, pf.lon_deg, pf.h_m, ell.a, ell.e2);
    const double load_s = wall_now() - t_load0;
    const std::string model_name = model.name.empty() ? a.gfc : model.name;
    std::printf("model %s  nmax=%d  tide_system=%s  points=%d (%s)  recurrence=%s  repeat=%d\n", model_name.c_str(),
                model.nmax, model.tide_system.c_str(), np, a.set.c_str(), a.recurrence.c_str(), a.repeat);
    std::printf("load_s=%.3f\n", load_s);

    if (a.threads > 0)
    {
#ifdef _OPENMP
      omp_set_num_threads(a.threads);
#endif
    }
    const int threads = a.threads > 0 ? a.threads : omp_threads_now();

    // Per-point modes: Legendre function and the recurrence tables (built once; shared by repetitions).
    LegendreFn legendre_fn = nullptr;
    if (a.recurrence == "sfc")
      legendre_fn = &legendre::sfc;
    else if (a.recurrence == "msfc")
      legendre_fn = &legendre::msfc;
    else if (a.recurrence == "era")
      legendre_fn = &legendre::era;
    Outputs last;
    for (int rep = 1; rep <= a.repeat; ++rep)
    {
      const double w0 = wall_now();
      const double c0 = cpu_now();

      Outputs o;
      if (a.recurrence == "fast")
      {
        o = run_fast(model, pts, a.nmax, a.threads);
      }
      else
      {
        // Timed: the disturbing model and the recurrence tables are part of the per-point mode, as the recurrence
        // tables are part of compute_functionals() in the fast mode.
        const Model dist = disturbing_model(model, ell);
        const legendre::Coefs coefs = legendre::make_coefs(a.nmax);
        o = run_per_point(model, dist, pts, a.nmax, ell, legendre_fn, coefs, a.threads);
      }

      const double wall = wall_now() - w0;
      const double cpu = cpu_now() - c0;
      std::printf("repeat %d/%d: wall=%.4f s  cpu=%.4f s\n", rep, a.repeat, wall, cpu);
      append_timing(a.timings, a, np, threads, rep, wall, cpu, load_s, model_name);
      last = std::move(o);
    }

    write_result(a.out, pf, last);
    write_companion_md(a.out, a, command, np, threads, load_s, model, model_name);
    std::printf("wrote %s (%d points) [V, T, dg]\n", a.out.c_str(), np);
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "gsh_shs_bench error: %s\n", e.what());
    return 1;
  }
}
