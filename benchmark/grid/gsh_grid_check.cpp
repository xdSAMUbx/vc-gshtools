// gsh_grid_check.cpp: accuracy check of the grid path against the point path at high degree (no external engine).
//
// For each case the grid is synthesised with gsh::synth_grid (the full grid, as in production); the selected rows are
// extracted and assembled into V, T and dg (grid_assembly.hpp). The same rows are then evaluated with the point path:
// grid_to_points() on the selected rows only, followed by gsh::compute_functionals() with outputs V, T and dg.
// Per field the check reports
//     max |grid - points| / max |points|        (criterion 1e-12, fixed beforehand; not to be changed)
// and the grid row of the largest absolute difference. The result is a Markdown table on stdout and in
// <out-dir>/grid_vs_points.md.
//
// Rows exactly at a pole (|lat| = 90 deg) are not selected: compute_functionals rejects dg at an exact pole.
// Rows near a pole (cos lat < 5e-3) are compared too: both paths use the full-triangle oracle there.
//
// Usage:
//   gsh_grid_check --gfc <file.gfc> [--out-dir D] [--cases gl,regional,dh1] [--nmax-high 2190] [--nmax-light 720]
//                  [--print-ref-quad]
// --print-ref-quad prints, without running them, the ref_quad command lines for the selected Gauss-Legendre rows
// (benchmark/meridian/ref_quad, built and run under WSL). ref_quad evaluates at longitude 0 only.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
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
  using Clock = std::chrono::steady_clock;

  constexpr double kCriterion = 1e-12;

  struct Args
  {
    std::string gfc, out_dir = ".";
    std::vector<std::string> cases = {"gl", "regional", "dh1"};
    int nmax_high = 2190, nmax_light = 720;
    bool print_ref_quad = false;
  };

  [[noreturn]] void usage(const char *msg)
  {
    if (msg)
      std::fprintf(stderr, "gsh_grid_check: %s\n", msg);
    std::fprintf(stderr, "usage: gsh_grid_check --gfc <file.gfc> [--out-dir D] [--cases gl,regional,dh1] [--nmax-high N] "
                         "[--nmax-light N] [--print-ref-quad]\n");
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
      else if (key == "--out-dir")
        a.out_dir = value();
      else if (key == "--cases")
      {
        a.cases.clear();
        const std::string s = value();
        std::size_t pos = 0;
        while (pos <= s.size())
        {
          const std::size_t comma = s.find(',', pos);
          a.cases.push_back(s.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos));
          if (comma == std::string::npos)
            break;
          pos = comma + 1;
        }
      }
      else if (key == "--nmax-high")
        a.nmax_high = std::atoi(value().c_str());
      else if (key == "--nmax-light")
        a.nmax_light = std::atoi(value().c_str());
      else if (key == "--print-ref-quad")
        a.print_ref_quad = true;
      else
        usage(("unknown option " + key).c_str());
    }
    if (a.gfc.empty())
      usage("--gfc is required");
    for (const auto &c : a.cases)
      if (c != "gl" && c != "regional" && c != "dh1" && c != "dh2")
        usage("--cases entries must be gl, regional, dh1 or dh2");
    return a;
  }

  struct FieldCmp
  {
    std::string name;
    double maxdiff = 0, maxabs = 0, rel = 0, lat_of_max = 0;
    std::int64_t row_of_max = -1;
    bool ok = false;
  };

  // Compares grid and point values point by point. Points are ordered as the selected rows, nlon per row.
  FieldCmp compare(const std::string &name, const std::vector<double> &g, const AlignedVec<double> &p, const std::vector<std::int64_t> &rows,
                   const Grid &grid)
  {
    FieldCmp c;
    c.name = name;
    const std::int64_t nlon = grid.nlon;
    const double inf = std::numeric_limits<double>::infinity();
    std::int64_t arg = 0;
    for (std::size_t k = 0; k < p.size(); ++k)
    {
      const double d = std::abs(g[k] - p[k]);
      if (std::isnan(d))
      {
        c.maxdiff = inf;
        arg = static_cast<std::int64_t>(k);
      }
      else if (d > c.maxdiff)
      {
        c.maxdiff = d;
        arg = static_cast<std::int64_t>(k);
      }
      if (!std::isnan(p[k]))
        c.maxabs = std::max(c.maxabs, std::abs(p[k]));
    }
    c.rel = c.maxdiff / c.maxabs;
    c.ok = c.rel <= kCriterion;
    c.row_of_max = rows[arg / nlon];
    c.lat_of_max = rad2deg(grid.lat[c.row_of_max]);
    return c;
  }

  std::string fmt(const char *f, double v)
  {
    char b[64];
    std::snprintf(b, sizeof b, f, v);
    return b;
  }

  struct CaseOut
  {
    std::string name;
    int nmax = 0;
    std::int64_t nrows = 0, npoints = 0;
    double grid_s = 0, points_s = 0;
    std::vector<FieldCmp> fields;
    std::vector<std::int64_t> gl_rows; // rows for the ref_quad command lines (GL only)
    std::vector<double> gl_lat_deg;
  };
}

int main(int argc, char **argv)
{
  try
  {
    const Args a = parse_args(argc, argv);
    const Ellipsoid ell = Ellipsoid::GRS80();
    const std::string gfc_name = a.gfc.substr(a.gfc.find_last_of("/\\") == std::string::npos ? 0 : a.gfc.find_last_of("/\\") + 1);
    const Model full = load_gfc(a.gfc, a.nmax_high);

    std::vector<CaseOut> outs;
    for (const std::string &name : a.cases)
    {
      CaseOut out;
      out.name = name;
      out.nmax = name == "dh1" || name == "dh2" ? a.nmax_light : a.nmax_high;
      const Model model = out.nmax == full.nmax ? full : full.truncated(out.nmax);
      const bool regional = name == "regional";
      const Grid grid = make_grid(name, out.nmax, model.R);

      // Rows: selection, minus rows at an exact pole.
      std::vector<std::int64_t> rows;
      for (std::int64_t i : parse_rows("auto", grid, regional))
        if (std::abs(grid.lat[i]) != PI / 2)
          rows.push_back(i);
      out.nrows = static_cast<std::int64_t>(rows.size());
      out.npoints = out.nrows * grid.nlon;

      // Grid side: full synthesis, then extraction of the selected rows.
      auto t0 = Clock::now();
      const Prepared pm = prepare(model, 0, out.nmax);
      const Model dist = disturbing_model(model, ell);
      const Prepared pd = prepare(dist, 0, out.nmax);
      RawFields rm = synth_grid(pm, grid, 0);
      RawFields rd = synth_grid(pd, grid, 1);
      std::vector<double> gV(out.npoints), gT(out.npoints), gDg(out.npoints);
      for (std::int64_t ri = 0; ri < out.nrows; ++ri)
      {
        const std::int64_t row = rows[ri];
        const double lat = grid.lat[row], r = grid.r[row];
        for (std::int64_t j = 0; j < grid.nlon; ++j)
        {
          const std::int64_t k = row * grid.nlon + j, p = ri * grid.nlon + j;
          gV[p] = potential(model.GM, r, rm.f[kS0][k]);
          gT[p] = potential(model.GM, r, rd.f[kS0][k]);
          gDg[p] = gravity_disturbance(ell, model.GM, lat, r, rd.f[kSR][k], rd.f[kSP][k], rd.f[kSL][k]);
        }
      }
      rm = RawFields{};
      rd = RawFields{};
      out.grid_s = std::chrono::duration<double>(Clock::now() - t0).count();

      // Point side: Points for the selected rows only, then the functional path.
      t0 = Clock::now();
      const Points pts = grid_to_points(sub_grid(grid, rows));
      FunctionalRequest req;
      req.outputs = {Functional::V, Functional::T, Functional::dg};
      req.nmax = out.nmax;
      const FunctionalResult res = compute_functionals(model, pts, req);
      out.points_s = std::chrono::duration<double>(Clock::now() - t0).count();

      out.fields.push_back(compare("V", gV, res.values.at(Functional::V), rows, grid));
      out.fields.push_back(compare("T", gT, res.values.at(Functional::T), rows, grid));
      out.fields.push_back(compare("dg", gDg, res.values.at(Functional::dg), rows, grid));

      if (name == "gl")
        for (std::int64_t row : rows)
        {
          out.gl_rows.push_back(row);
          out.gl_lat_deg.push_back(rad2deg(grid.lat[row]));
        }
      std::fprintf(stderr, "case %s: nmax=%d rows=%lld points=%lld grid=%.1f s points=%.1f s\n", name.c_str(), out.nmax,
                   static_cast<long long>(out.nrows), static_cast<long long>(out.npoints), out.grid_s, out.points_s);
      outs.push_back(std::move(out));
    }

    // ---- Markdown report ----------------------------------------------------------------------------------------
    std::string md;
    md += "# Grid path vs point path (gsh)\n\n";
    md += "Criterion: max|grid - points| / max|points| <= 1e-12 for each field (fixed before the run).\n\n";
    md += "Model: `" + gfc_name + "`, truncated to the degree of each case. Spherical latitudes, constant radius R per row.\n";
    md += "Grid side: `synth_grid` on the full grid (V: grad 0 of the model; T and dg: grad 1 of the disturbing model), assembled\n";
    md += "with the formulas of `functionals.cpp`. Point side: `grid_to_points` on the selected rows, then `compute_functionals`\n";
    md += "with outputs V, T, dg. Selected rows: auto (global) or the 5 regional rows; exact-pole rows excluded.\n\n";
    md += "| case | nmax | rows | points | field | max abs diff | max abs value | relative diff | row of max diff (lat deg) | verdict |\n";
    md += "|---|---|---|---|---|---|---|---|---|---|\n";
    for (const auto &o : outs)
      for (const auto &f : o.fields)
        md += "| " + o.name + " | " + std::to_string(o.nmax) + " | " + std::to_string(o.nrows) + " | " + std::to_string(o.npoints) +
              " | " + f.name + " | " + fmt("%.3e", f.maxdiff) + " | " + fmt("%.3e", f.maxabs) + " | " + fmt("%.3e", f.rel) +
              " | row " + std::to_string(f.row_of_max) + " (" + fmt("%.6f", f.lat_of_max) + ") | " + (f.ok ? "PASS" : "FAIL") + " |\n";
    md += "\nWall time of each side in this single accuracy run (not a timing benchmark):\n\n";
    md += "| case | grid side (s) | point side (s) |\n|---|---|---|\n";
    for (const auto &o : outs)
      md += "| " + o.name + " | " + fmt("%.1f", o.grid_s) + " | " + fmt("%.1f", o.points_s) + " |\n";

    if (a.print_ref_quad)
    {
      for (const auto &o : outs)
        if (o.name == "gl" && !o.gl_lat_deg.empty())
        {
          std::string list;
          for (std::size_t i = 0; i < o.gl_lat_deg.size(); ++i)
          {
            if (i)
              list += ",";
            list += fmt("%.17g", o.gl_lat_deg[i]);
          }
          md += "\nref_quad command line for the selected GL rows (NOT run; WSL; engine spherical latitude in degrees):\n\n";
          md += "```\nref_quad " + gfc_name + " ref_quad_gl_rows.csv --nmax " + std::to_string(o.nmax) + " --engine-rad --lat-list " + list +
                "\n```\n";
        }
    }

    std::fputs(md.c_str(), stdout);
    std::filesystem::create_directories(a.out_dir);
    const std::string path = a.out_dir + "/grid_vs_points.md";
    std::FILE *f = std::fopen(path.c_str(), "w");
    if (!f)
      throw std::runtime_error("cannot write " + path);
    std::fputs(md.c_str(), f);
    std::fclose(f);
    std::fprintf(stderr, "wrote %s\n", path.c_str());
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "gsh_grid_check: %s\n", e.what());
    return 1;
  }
}
