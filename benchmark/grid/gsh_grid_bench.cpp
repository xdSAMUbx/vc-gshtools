// gsh_grid_bench.cpp: gsh driver of the grid benchmark (definitions: benchmark/grid/README.md, section "gsh side").
//
// One run = one case (gl, dh1, dh2 or regional) at one degree. The grid is synthesised K times with gsh::synth_grid:
//   V  from the model (grad 0, field S_0);
//   T  from the disturbing model (grad 1 when dg is requested, which also gives S_0; grad 0 otherwise);
//   dg from the disturbing gradient sums (grad 1) and the normal gravity of the row (Ellipsoid::normal_gravity).
// Each repetition appends one line to <out-dir>/timings.csv. The functionals are assembled in memory for the selected
// rows only; the full grid is never written to disk. The selected rows go to <out-dir>/gsh_<case>_rows.csv.
//
// Usage:
//   gsh_grid_bench --gfc <file.gfc> --nmax N --case gl|dh1|dh2|regional --recurrence fast|msfc --threads T
//                  --repeat K --out-dir D --rows auto|i,j,k --fields V,T,dg
//
// Timings: wall_s is the sum of GridTiming.total over the synth_grid calls of one repetition (the assembly of the
// selected rows is not timed); lumps_thread_s and fft_thread_s are the matching sums of the thread-seconds.
// peak_mem_mb is the process peak working set (GetProcessMemoryInfo, Windows), read after each repetition, so it
// is cumulative over the process.
//
// --recurrence msfc is rejected: the grid path (synth_grid) has only the fast kernel tables. The MSFC recurrence exists
// only in the per-point oracle (reference_point), which is too slow for grids of this size.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif
#include "ellipsoid.hpp"
#include "functionals.hpp"
#include "grid.hpp"
#include "model.hpp"
#include "shs.hpp"
#include "grid_assembly.hpp"

namespace
{
  using namespace gsh;
  using namespace gridbench;

  struct Args
  {
    std::string gfc, case_name, recurrence = "fast", out_dir = ".", rows = "auto", fields = "V,T,dg";
    int nmax = -1, threads = 0, repeat = 1;
  };

  [[noreturn]] void usage(const char *msg)
  {
    if (msg)
      std::fprintf(stderr, "gsh_grid_bench: %s\n", msg);
    std::fprintf(stderr,
                 "usage: gsh_grid_bench --gfc <file.gfc> --nmax N --case gl|dh1|dh2|regional --recurrence fast|msfc\n"
                 "                      --threads T --repeat K --out-dir D --rows auto|i,j,k --fields V,T,dg\n");
    std::exit(2);
  }

  Args parse_args(int argc, char **argv)
  {
    Args a;
    for (int i = 1; i < argc; ++i)
    {
      const std::string key = argv[i];
      const auto value = [&]() -> std::string
      {
        if (i + 1 >= argc)
          usage(("missing value after " + key).c_str());
        return argv[++i];
      };
      if (key == "--gfc")
        a.gfc = value();
      else if (key == "--nmax")
        a.nmax = std::atoi(value().c_str());
      else if (key == "--case")
        a.case_name = value();
      else if (key == "--recurrence")
        a.recurrence = value();
      else if (key == "--threads")
        a.threads = std::atoi(value().c_str());
      else if (key == "--repeat")
        a.repeat = std::atoi(value().c_str());
      else if (key == "--out-dir")
        a.out_dir = value();
      else if (key == "--rows")
        a.rows = value();
      else if (key == "--fields")
        a.fields = value();
      else
        usage(("unknown option " + key).c_str());
    }
    if (a.gfc.empty() || a.case_name.empty() || a.nmax < 0)
      usage("--gfc, --case and --nmax are required");
    if (a.recurrence != "fast" && a.recurrence != "msfc")
      usage("--recurrence must be fast or msfc");
    if (a.repeat < 1 || a.threads < 0)
      usage("--repeat must be >= 1 and --threads >= 0");
    return a;
  }

  bool has_field(const std::string &fields, const char *name)
  {
    std::size_t pos = 0;
    while (pos <= fields.size())
    {
      const std::size_t comma = fields.find(',', pos);
      if (fields.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos) == name)
        return true;
      if (comma == std::string::npos)
        break;
      pos = comma + 1;
    }
    return false;
  }

  // Process peak working set in MiB (Windows). NaN elsewhere.
  double peak_mem_mb()
  {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc{};
    pmc.cb = sizeof pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc))
      return static_cast<double>(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0);
#endif
    return std::numeric_limits<double>::quiet_NaN();
  }

  void append_timing(const std::string &path, const Args &a, const std::string &recurrence_label, int threads, std::int64_t nlat,
                     std::int64_t nlon, int rep, double wall, double lumps, double fft, double mem)
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
      f << "engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,lumps_thread_s,fft_thread_s,peak_mem_mb\n";
    char buf[512];
    std::snprintf(buf, sizeof buf, "gsh,%s,%s,%d,%d,%lld,%lld,%d,%.6f,%.6f,%.6f,%.1f\n", a.case_name.c_str(), recurrence_label.c_str(), threads,
                  a.nmax, static_cast<long long>(nlat), static_cast<long long>(nlon), rep, wall, lumps, fft, mem);
    f << buf;
  }
}

int main(int argc, char **argv)
{
  try
  {
    const Args a = parse_args(argc, argv);
    if (a.recurrence == "msfc")
      throw std::invalid_argument("--recurrence msfc is not available on the grid path: synth_grid uses the fast kernel "
                                  "tables only (MSFC exists in the per-point oracle, too slow for grids)");
    const bool want_V = has_field(a.fields, "V"), want_T = has_field(a.fields, "T"), want_dg = has_field(a.fields, "dg");
    if (!want_V && !want_T && !want_dg)
      throw std::invalid_argument("--fields must name at least one of V, T, dg");
    const bool regional = a.case_name == "regional";
    const Ellipsoid ell = Ellipsoid::GRS80();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // ---- Load, grid and selection (not timed) ----------------------------------------------------------------
    const Model model = load_gfc(a.gfc, a.nmax);
    const Grid grid = make_grid(a.case_name, a.nmax, model.R);
    const std::vector<std::int64_t> rows = parse_rows(a.rows, grid, regional);
    const Prepared prep_model = prepare(model, 0, a.nmax);
    Model dist;
    Prepared prep_dist;
    if (want_T || want_dg)
    {
      dist = disturbing_model(model, ell);
      prep_dist = prepare(dist, 0, a.nmax);
    }
    const int grad_dist = want_dg ? 1 : 0;
    const int threads = a.threads > 0 ? a.threads : omp_get_max_threads();
    SynthOptions opts;
    opts.threads = a.threads;
    const std::string timings = a.out_dir + "/timings.csv";

    // ---- Repetitions -----------------------------------------------------------------------------------------
    RawFields raw_model, raw_dist;
    for (int rep = 1; rep <= a.repeat; ++rep)
    {
      double wall = 0, lumps = 0, fft = 0;
      if (want_V)
      {
        raw_model = RawFields{};
        GridTiming tm;
        raw_model = synth_grid(prep_model, grid, 0, opts, &tm);
        wall += tm.total;
        lumps += tm.lumps;
        fft += tm.fft;
      }
      if (want_T || want_dg)
      {
        raw_dist = RawFields{};
        GridTiming td;
        raw_dist = synth_grid(prep_dist, grid, grad_dist, opts, &td);
        wall += td.total;
        lumps += td.lumps;
        fft += td.fft;
      }
      const double mem = peak_mem_mb();
      append_timing(timings, a, a.recurrence, threads, grid.nlat(), grid.nlon, rep, wall, lumps, fft, mem);
      std::printf("repeat %d/%d: wall=%.4f s  lumps_thread=%.4f s  fft_thread=%.4f s  peak=%.1f MiB\n", rep, a.repeat, wall, lumps,
                  fft, mem);
    }

    // ---- Assembly of the selected rows (last repetition) and output ------------------------------------------
    const std::string rows_path = a.out_dir + "/gsh_" + a.case_name + "_rows.csv";
    std::FILE *f = std::fopen(rows_path.c_str(), "w");
    if (!f)
      throw std::runtime_error("cannot write " + rows_path);
    std::fprintf(f, "row,lat_deg,lon_deg,V_m2s2,T_m2s2,dg_ms2\n");
    const double threshold = FunctionalRequest{}.pole_cos_threshold;
    for (std::int64_t row : rows)
    {
      const double lat = grid.lat[row], r = grid.r[row];
      const bool dg_ok = std::cos(lat) >= threshold;
      for (std::int64_t j = 0; j < grid.nlon; ++j)
      {
        const std::int64_t k = row * grid.nlon + j;
        const double lon = grid.lon0 + static_cast<double>(j) * grid.dlon;
        const double V = want_V ? potential(model.GM, r, raw_model.f[kS0][k]) : nan;
        const double T = want_T ? potential(model.GM, r, raw_dist.f[kS0][k]) : nan;
        const double dg = (want_dg && dg_ok) ? gravity_disturbance(ell, model.GM, lat, r, raw_dist.f[kSR][k], raw_dist.f[kSP][k], raw_dist.f[kSL][k]) : nan;
        std::fprintf(f, "%lld,%.17g,%.17g,%.17g,%.17g,%.17g\n", static_cast<long long>(row), rad2deg(lat), rad2deg(lon), V, T, dg);
      }
    }
    std::fclose(f);
    std::printf("wrote %s (%zu rows) and appended %s\n", rows_path.c_str(), rows.size(), timings.c_str());
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "gsh_grid_bench: %s\n", e.what());
    return 1;
  }
}
