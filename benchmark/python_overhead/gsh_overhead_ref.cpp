// gsh_overhead_ref.cpp: C++ reference of the Python-overhead benchmark (benchmark/python_overhead/README.md).
//
// For each job (a points file and a set of outputs) the driver times the same work that one Python call
// gshtools.compute_functionals() does, with the same inputs:
//   points  : gsh::Points::from_geodetic(lat, lon, h) from arrays already in memory
//   compute : gsh::compute_functionals(model, points, request)
// The model is loaded once, and the points file is read before the timed block; neither is timed.
// Wall time is std::chrono::steady_clock. The request (outputs, nmin, nmax, threads) is built before the timed block,
// as on the Python side.
//
// Usage:
//   gsh_overhead_ref --gfc <file.gfc> --nmax N --threads K --warmup W --repeat R
//                    --job <V|VTdg> <points.csv> [--job <V|VTdg> <points.csv> ...]
//                    [--dump-dir <dir>]
//   threads = 0 means the OpenMP default. --warmup W runs W untimed-by-the-statistics calls first (phase "warmup");
//   --repeat R timed calls (phase "timed"). --dump-dir writes, per job, the output arrays of the LAST call as raw
//   little-endian float64 (outputs in request order, one array after the other), for the bitwise check.
//
// Output (stdout), one line per call:
//   RUN,<job index>,<phase>,<rep>,<points_s>,<compute_s>,<stage_sum_s>
//   stage_sum_s is the sum of FunctionalResult.timing_s (the stages inside the library).
//   Other lines: CONFIG,... (threads, load time). Errors go to stderr with exit code 1.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "ellipsoid.hpp"
#include "functionals.hpp"
#include "model.hpp"
#include "synth_types.hpp"
#include "table_io.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr,
                 "gsh_overhead_ref: %s\n"
                 "usage: gsh_overhead_ref --gfc <file.gfc> --nmax N --threads K --warmup W --repeat R\n"
                 "                        --job <V|VTdg> <points.csv> [--job ...] [--dump-dir <dir>]\n",
                 why);
    std::exit(2);
  }

  struct Job
  {
    std::string name;   // "V" or "VTdg"
    std::string points; // CSV with header id,lat_deg,lon_deg,h_m
  };

  struct Args
  {
    std::string gfc, dump_dir;
    int nmax = -1, threads = 0, warmup = 1, repeat = 8;
    std::vector<Job> jobs;
  };

  int parse_int(const std::string &s, const char *what)
  {
    char *end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0')
      usage(what);
    return static_cast<int>(v);
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
      else if (key == "--nmax")
        a.nmax = parse_int(value(), "--nmax must be an integer");
      else if (key == "--threads")
        a.threads = parse_int(value(), "--threads must be an integer");
      else if (key == "--warmup")
        a.warmup = parse_int(value(), "--warmup must be an integer");
      else if (key == "--repeat")
        a.repeat = parse_int(value(), "--repeat must be an integer");
      else if (key == "--dump-dir")
        a.dump_dir = value();
      else if (key == "--job")
      {
        if (i + 2 >= argc)
          usage("--job needs two values: <V|VTdg> <points.csv>");
        Job j;
        j.name = argv[i + 1];
        j.points = argv[i + 2];
        i += 2;
        if (j.name != "V" && j.name != "VTdg")
          usage("--job functionals must be V or VTdg");
        a.jobs.push_back(j);
      }
      else
        usage(("unknown option " + key).c_str());
    }
    if (a.gfc.empty() || a.jobs.empty())
      usage("--gfc and at least one --job are required");
    if (a.nmax < 2)
      usage("--nmax must be at least 2");
    if (a.threads < 0 || a.warmup < 0 || a.repeat < 0)
      usage("--threads, --warmup and --repeat must not be negative");
    return a;
  }

  std::vector<gsh::Functional> outputs_of(const std::string &name)
  {
    if (name == "V")
      return {gsh::Functional::V};
    return {gsh::Functional::V, gsh::Functional::T, gsh::Functional::dg};
  }

  double wall_now()
  {
    using Clock = std::chrono::steady_clock;
    static const Clock::time_point origin = Clock::now();
    return std::chrono::duration<double>(Clock::now() - origin).count();
  }

  int omp_threads_now()
  {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    const Args a = parse_args(argc, argv);
    const gsh::Ellipsoid ell = gsh::Ellipsoid::GRS80();

    const double t_load0 = wall_now();
    const gsh::Model model = gsh::load_gfc(a.gfc, a.nmax);
    const double load_s = wall_now() - t_load0;
    std::printf("CONFIG,nmax=%d,threads_requested=%d,threads_openmp_max=%d,load_s=%.6f\n", a.nmax, a.threads,
                omp_threads_now(), load_s);

    for (std::size_t j = 0; j < a.jobs.size(); ++j)
    {
      const Job &job = a.jobs[j];

      // Points file: read before the timed block (same parser as gshtools.read_table, columns 1..3).
      const gsh::Table tab = gsh::read_table(job.points, ',', 1, {1, 2, 3});
      if (!tab.skipped_lines.empty())
        throw std::runtime_error("non-numeric line in " + job.points);
      const std::vector<double> &lat = tab.columns[0], &lon = tab.columns[1], &h = tab.columns[2];

      gsh::FunctionalRequest rq;
      rq.outputs = outputs_of(job.name);
      rq.nmin = 0;
      rq.nmax = a.nmax;
      rq.ellipsoid = ell;
      rq.synth.threads = a.threads;

      const int total_calls = a.warmup + a.repeat;
      gsh::FunctionalResult last;
      for (int call = 0; call < total_calls; ++call)
      {
        const bool warm = call < a.warmup;
        const int rep = warm ? 0 : call - a.warmup + 1;

        const double t0 = wall_now();
        const gsh::Points pts = gsh::Points::from_geodetic(lat, lon, h, ell.a, ell.e2);
        const double t1 = wall_now();
        gsh::FunctionalResult res = gsh::compute_functionals(model, pts, rq);
        const double t2 = wall_now();

        double stage_sum = 0.0;
        for (const auto &kv : res.timing_s)
          stage_sum += kv.second;
        std::printf("RUN,%zu,%s,%d,%.9f,%.9f,%.9f\n", j, warm ? "warmup" : "timed", rep, t1 - t0, t2 - t1, stage_sum);
        std::fflush(stdout);
        last = std::move(res);
      }

      if (!a.dump_dir.empty() && total_calls > 0)
      {
        const std::string path = a.dump_dir + "/job" + std::to_string(j) + "_threads" + std::to_string(a.threads) +
                                 "_" + job.name + ".bin";
        std::ofstream out(path, std::ios::binary);
        if (!out)
          throw std::runtime_error("cannot write " + path);
        for (gsh::Functional f : rq.outputs)
        {
          const auto &v = last.values.at(f);
          out.write(reinterpret_cast<const char *>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(double)));
        }
      }
    }
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "gsh_overhead_ref: %s\n", e.what());
    return 1;
  }
}
