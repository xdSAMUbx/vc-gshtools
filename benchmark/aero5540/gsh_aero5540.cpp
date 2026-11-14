// gsh_aero5540: GRAV-D flight (Colorado), all 45 functionals at degree/order 5540 with XGM2019e.
//
//   gsh_aero5540 <model.gfc> <flight.txt> <dem.grd> <out_dir> [nmax=5540] [reps=8] [sub_points=10000] [acc_points=200]
//
// Three timing measurements (each repeated `reps` times; the first repetition is NOT discarded, all are recorded):
//   full  : the 45 functionals together on every flight point. The library plans the minimum number of synthesis
//           passes, so this is the real cost of "everything at once". Stage times come from FunctionalResult::timing_s.
//   alone : each functional requested by itself on a fixed subsample (every k-th point, `sub_points` points). The time
//           of a functional alone includes the synthesis pass it needs; functionals of one group repeat the same pass,
//           so the per-functional times are NOT additive (the full run is cheaper than their sum).
//   sub45 : the 45 together on the same subsample, used to relate "alone" times to the full run.
// Accuracy (one pass, not repeated): V and the six local tensor components V_xx..V_zz at `acc_points` flight points,
// compared with the sequential oracle (reference_point, double precision) at the same degree. Metric per quantity:
// max|d| / max|ref| (the criterion of benchmark/meridian is 1e-12).
//
// Flight data: GRAV-D (NGS); coordinates only are used, the data are not redistributed. N and zeta use the Colorado
// DEM (5' mean heights); points outside the DEM give NaN for those two functionals.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif
#include "functionals.hpp"
#include "table_io.hpp"
#include "shs.hpp"
#include "legendre.hpp"
using namespace gsh;
using F = Functional;

namespace
{
  double secs(std::chrono::steady_clock::time_point a) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count(); }

  double peak_mem_mb()
  {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc))
      return static_cast<double>(pmc.PeakWorkingSetSize) / 1048576.0;
#endif
    return -1.0;
  }

  Points subsample(const Points &p, std::int64_t k)
  {
    Points s;
    for (std::int64_t i = 0; i < p.size(); i += k)
    {
      s.lat.push_back(p.lat[i]);
      s.lon.push_back(p.lon[i]);
      s.r.push_back(p.r[i]);
    }
    return s;
  }
}

int main(int argc, char **argv)
{
  if (argc < 5)
  {
    std::fprintf(stderr, "usage: %s <model.gfc> <flight.txt> <dem.grd> <out_dir> [nmax=5540] [reps=8] [sub_points=10000] [acc_points=200]\n", argv[0]);
    return 2;
  }
  const int nmax = argc > 5 ? std::atoi(argv[5]) : 5540;
  const int reps = argc > 6 ? std::atoi(argv[6]) : 8;
  const std::int64_t sub_n = argc > 7 ? std::atoll(argv[7]) : 10000;
  const std::int64_t acc_n = argc > 8 ? std::atoll(argv[8]) : 200;
  // Optional: repetitions of the all-points run (default = reps) and how many points it uses (0 = all of them).
  const int full_reps = argc > 9 ? std::atoi(argv[9]) : reps;
  const std::int64_t full_n = argc > 10 ? std::atoll(argv[10]) : 0;
  const bool resume = argc > 11 && std::atoi(argv[11]) != 0;
  // Optional (12): 1 = fused-pass A/B: the all-points run is repeated with fuse_model_and_disturbing true and false and the two
  // results are compared. Needs full_reps > 0. sub_points = 0 skips the per-functional stage, acc_points = 0 the oracle stage.
  const bool ab = argc > 12 && std::atoi(argv[12]) != 0;
  const std::string out = argv[4];
  std::setvbuf(stdout, nullptr, _IONBF, 0); // progress must be visible when stdout is redirected to a file
  std::filesystem::create_directories(out);
  const Ellipsoid ell = Ellipsoid::GRS80();

  auto t0 = std::chrono::steady_clock::now();
  const Model model = load_gfc(argv[1], nmax);
  std::printf("model nmax=%d read in %.1f s (%s)\n", nmax, secs(t0), model.tide_system.c_str());
  const auto fl = read_table(argv[2], ',', 1, {3, 4, 5});
  const Points all_pts = Points::from_geodetic(fl.columns[0], fl.columns[1], fl.columns[2], ell.a, ell.e2);
  const Points pts = (full_n > 0 && full_n < all_pts.size()) ? subsample(all_pts, std::max<std::int64_t>(1, all_pts.size() / full_n)) : all_pts;
  std::printf("flight: %lld points\n", static_cast<long long>(pts.size()));
  auto dem = load_gravsoft_grid(argv[3]);

  FunctionalRequest rq;
  rq.outputs = all_functionals();
  rq.nmax = nmax;
  rq.height = height_from_grid(std::move(dem), 5.0 / 60.0);
  std::printf("ISA: %s\n", active_isa_name().c_str());

  std::FILE *ft = full_reps > 0 ? std::fopen((out + "/timings_full.csv").c_str(), "w") : nullptr; // no empty file when the all-together run is skipped
  if (ft)
    std::fprintf(ft, "fuse,rep,points,total_s,model_synthesis_s,disturbing_synthesis_s,height_s,ellipsoid_synthesis_s,polar_points,points_without_height,peak_mem_mb\n");
  FunctionalResult last, last_unfused;
  for (int mode = 0; mode < (ab ? 2 : 1); ++mode)
  {
    rq.fuse_model_and_disturbing = (mode == 0);
    for (int r = 0; r < full_reps; ++r)
    {
      t0 = std::chrono::steady_clock::now();
      FunctionalResult res = compute_functionals(model, pts, rq);
      const double tot = secs(t0);
      auto g = [&](const char *k) { auto it = res.timing_s.find(k); return it == res.timing_s.end() ? 0.0 : it->second; };
      std::fprintf(ft, "%d,%d,%lld,%.4f,%.4f,%.4f,%.4f,%.4f,%lld,%lld,%.1f\n", rq.fuse_model_and_disturbing ? 1 : 0, r, static_cast<long long>(pts.size()), tot, g("model_synthesis"), g("disturbing_synthesis"), g("height"), g("ellipsoid_synthesis"), static_cast<long long>(res.polar_points), static_cast<long long>(res.points_without_height), peak_mem_mb());
      std::fflush(ft);
      std::printf("full fuse=%d rep %d: %.1f s\n", rq.fuse_model_and_disturbing ? 1 : 0, r, tot);
      if (r + 1 == full_reps)
        (mode == 0 ? last : last_unfused) = std::move(res);
    }
  }
  rq.fuse_model_and_disturbing = true;
  if (ft)
    std::fclose(ft);
  if (ab && full_reps > 0)
  {
    // fused vs unfused, same points: max|d| / max|ref| per functional (ref = unfused two-pass result)
    std::FILE *fd = std::fopen((out + "/fuse_vs_unfused.csv").c_str(), "w");
    std::fprintf(fd, "functional,max_abs_diff,max_abs_ref,relative\n");
    double worst = 0;
    for (F f : all_functionals())
    {
      const auto &a = last.values.at(f);
      const auto &b = last_unfused.values.at(f);
      double md = 0, mr = 0;
      for (std::size_t i = 0; i < a.size(); ++i)
        if (std::isfinite(a[i]) && std::isfinite(b[i]))
        {
          md = std::max(md, std::fabs(a[i] - b[i]));
          mr = std::max(mr, std::fabs(b[i]));
        }
      std::fprintf(fd, "%s,%.6e,%.6e,%.6e\n", functional_name(f), md, mr, mr > 0 ? md / mr : 0.0);
      worst = std::max(worst, mr > 0 ? md / mr : 0.0);
    }
    std::fclose(fd);
    std::printf("fused vs unfused: worst max|d|/max|ref| = %.3e\n", worst);
  }

  // ---- each functional alone on the subsample ----
  const std::int64_t k = std::max<std::int64_t>(1, pts.size() / std::max<std::int64_t>(1, sub_n));
  const Points sub = subsample(pts, k);
  // resume (argument 11 = 1): keep the rows already in timings_alone.csv, skip those functionals and append the rest.
  std::set<std::string> already_done;
  const std::string alone_path = out + "/timings_alone.csv";
  if (resume)
  {
    std::ifstream in(alone_path);
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
      if (!line.empty())
        already_done.insert(line.substr(0, line.find(',')));
    std::printf("resume: %zu functionals already in %s\n", already_done.size(), alone_path.c_str());
  }
  std::FILE *fa = std::fopen(alone_path.c_str(), resume ? "a" : "w");
  if (!resume)
    std::fprintf(fa, "functional,rep,points,wall_s\n");
  for (int r = 0; r < (full_reps > 0 ? reps : 0); ++r) // ALL45 on the subsample only when the all-together run is wanted
  {
    FunctionalRequest all = rq;
    t0 = std::chrono::steady_clock::now();
    (void)compute_functionals(model, sub, all);
    std::fprintf(fa, "ALL45,%d,%lld,%.4f\n", r, static_cast<long long>(sub.size()), secs(t0));
  }
  for (F f : (sub_n > 0 ? all_functionals() : std::vector<F>{}))
  {
    if (already_done.count(functional_name(f)))
      continue;
    FunctionalRequest one = rq;
    one.outputs = {f};
    for (int r = 0; r < reps; ++r)
    {
      t0 = std::chrono::steady_clock::now();
      (void)compute_functionals(model, sub, one);
      const double t = secs(t0);
      std::fprintf(fa, "%s,%d,%lld,%.4f\n", functional_name(f), r, static_cast<long long>(sub.size()), t);
      std::fflush(fa);
    }
    std::printf("alone %-16s done\n", functional_name(f));
  }
  std::fclose(fa);

  if (acc_n <= 0)
    return 0;
  // ---- accuracy against the oracle (single pass) ----
  const legendre::Coefs coefs = legendre::make_coefs(nmax);
  const std::int64_t ka = std::max<std::int64_t>(1, pts.size() / std::max<std::int64_t>(1, acc_n));
  struct Acc { const char *name; F f; double maxd = 0, maxr = 0; };
  std::vector<Acc> acc = {{"V", F::V}, {"V_xx", F::V_xx}, {"V_xy", F::V_xy}, {"V_xz", F::V_xz}, {"V_yy", F::V_yy}, {"V_yz", F::V_yz}, {"V_zz", F::V_zz}};
  std::int64_t used = 0;
  // The accuracy points are computed by themselves (7 outputs), so this stage does not need the all-points run.
  const Points accpts = subsample(pts, ka);
  FunctionalRequest racc = rq;
  racc.outputs = {F::V, F::V_xx, F::V_xy, F::V_xz, F::V_yy, F::V_yz, F::V_zz};
  const FunctionalResult accres = compute_functionals(model, accpts, racc);
  for (std::int64_t i = 0; i < accpts.size(); ++i)
  {
    const ReferencePoint ref = reference_point(model, coefs, accpts.lat[i], accpts.lon[i], accpts.r[i], 0, nmax);
    const double rv[7] = {model.GM / accpts.r[i] * ref.raw[0], ref.vxx, ref.vxy, ref.vxz, ref.vyy, ref.vyz, ref.vzz};
    for (std::size_t j = 0; j < acc.size(); ++j)
    {
      const double v = accres.values.at(acc[j].f)[i];
      acc[j].maxd = std::max(acc[j].maxd, std::fabs(v - rv[j]));
      acc[j].maxr = std::max(acc[j].maxr, std::fabs(rv[j]));
    }
    ++used;
  }
  std::FILE *fc = std::fopen((out + "/accuracy.csv").c_str(), "w");
  std::fprintf(fc, "functional,points,max_abs_diff,max_abs_ref,relative\n");
  for (const auto &a : acc)
  {
    std::fprintf(fc, "%s,%lld,%.6e,%.6e,%.6e\n", a.name, static_cast<long long>(used), a.maxd, a.maxr, a.maxd / a.maxr);
    std::printf("accuracy %-5s max|d|/max|ref| = %.3e\n", a.name, a.maxd / a.maxr);
  }
  std::fclose(fc);
  return 0;
}
