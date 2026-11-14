// shs.cpp: scheduler for scattered-point synthesis. It selects the instruction set (kernels::Table), the tile
// width and the parallel strategy, then runs the kernel on every tile of points. The arithmetic lives in
// gsh_legendre/kernels; this file only prepares inputs and gathers outputs.
#include "shs.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <utility>
#include "kernels.hpp"
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace gsh
{
  // Builds the recurrence tables once; a and b are interleaved so that the kernel reads one contiguous pair per (n, m).
  std::shared_ptr<const Recurrence> make_recurrence(int nmax)
  {
    auto rec = std::make_shared<Recurrence>();
    rec->nmax = nmax;
    auto c = legendre::make_coefs(nmax);
    const auto size = tri_size(nmax);
    rec->ab.resize(2 * size);
    // Moved, not copied: c is local and only a and b are read below. Same values, no serial copy of the triangles.
    rec->f = std::move(c.f);
    rec->d = std::move(c.d);
#pragma omp parallel for if (size > 100000)
    for (std::int64_t i = 0; i < size; ++i)
    {
      rec->ab[2 * i] = c.a[i];
      rec->ab[2 * i + 1] = c.b[i];
    }
    return rec;
  }

  // Repacks C and S from the model's triangle (model.nmax) into the triangle of nmax, zeroing n < nmin.
  // The copy is parallel over m; the packed layout is the one described in legendre.hpp.
  Prepared prepare(const Model &model, int nmin, int nmax, std::shared_ptr<const Recurrence> rec)
  {
    if (nmax < 0)
      nmax = model.nmax;
    if (nmax > model.nmax)
      throw std::invalid_argument("prepare: nmax is larger than the model degree");
    if (nmin < 0 || nmin > nmax)
      throw std::invalid_argument("prepare: nmin is out of range");
    if (rec && rec->nmax != nmax)
      throw std::invalid_argument("prepare: the recurrence tables are for another degree");
    Prepared p;
    p.nmin = nmin;
    p.nmax = nmax;
    p.GM = model.GM;
    p.R = model.R;
    p.rec = rec ? std::move(rec) : make_recurrence(nmax);
    p.cs.resize(2 * tri_size(nmax));
#pragma omp parallel for schedule(dynamic, 16) if (tri_size(nmax) > 100000)
    for (int m = 0; m <= nmax; ++m)
    {
      const auto dst = col_offset(m, nmax), src = col_offset(m, model.nmax);
      for (int n = m; n <= nmax; ++n)
      {
        const auto i = dst + (n - m), j = src + (n - m);
        p.cs[2 * i] = n >= nmin ? model.C[j] : 0.0;
        p.cs[2 * i + 1] = n >= nmin ? model.S[j] : 0.0;
      }
    }
    return p;
  }

  namespace
  {
    // AVX-512 is usable only if the CPU reports AVX512F, DQ, VL and BW (CPUID leaf 7) and the OS saves the
    // opmask and ZMM state (XCR0 bits 1, 2, 5, 6, 7 = mask 0xE6, read with XGETBV). The list matches the instruction
    // sets the AVX-512 kernel is compiled with (/arch:AVX512 or -mavx512f/dq/vl/bw), so a CPU that passes this check
    // can run every instruction in kernels_avx512.cpp.
    bool cpu_has_avx512()
    {
#if defined(_MSC_VER)
      int info[4];
      __cpuidex(info, 7, 0);
      const bool f = (info[1] >> 16) & 1, dq = (info[1] >> 17) & 1, bw = (info[1] >> 30) & 1, vl = (info[1] >> 31) & 1;
      __cpuid(info, 1);
      const bool osxsave = (info[2] >> 27) & 1;
      if (!osxsave)
        return false;
      const unsigned long long xcr = _xgetbv(0);
      return f && dq && vl && bw && (xcr & 0xE6) == 0xE6;
#else
      return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
             __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw");
#endif
    }
    // AVX2 (CPUID leaf 7, EBX bit 5) and FMA (CPUID leaf 1, ECX bit 12), plus OS support for YMM state (XCR0 mask 6).
    bool cpu_has_avx2()
    {
#if defined(_MSC_VER)
      int info[4];
      __cpuidex(info, 7, 0);
      const bool avx2 = (info[1] >> 5) & 1;
      __cpuid(info, 1);
      const bool fma = (info[2] >> 12) & 1, osxsave = (info[2] >> 27) & 1;
      if (!osxsave)
        return false;
      return avx2 && fma && (_xgetbv(0) & 6) == 6;
#else
      return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#endif
    }
    // Kernel table for the requested ISA. An explicit request for an unsupported ISA throws; Auto picks AVX-512,
    // then AVX2, then the scalar table.
    const kernels::Table &select(Isa isa)
    {
      switch (isa)
      {
      case Isa::Scalar:
        return kernels::table_scalar();
      case Isa::Avx2:
        if (!cpu_has_avx2())
          throw std::runtime_error("CPU has no AVX2/FMA");
        return kernels::table_avx2();
      case Isa::Avx512:
        if (!cpu_has_avx512())
          throw std::runtime_error("CPU has no AVX-512");
        return kernels::table_avx512();
      default:
        if (cpu_has_avx512())
          return kernels::table_avx512();
        if (cpu_has_avx2())
          return kernels::table_avx2();
        return kernels::table_scalar();
      }
    }
    // Maps the number of SIMD vectors per tile (1, 2, 4, 8, 16) to the table index 0..4.
    int b_index(int vectors)
    {
      switch (vectors)
      {
      case 1:
        return 0;
      case 2:
        return 1;
      case 4:
        return 2;
      case 8:
        return 3;
      case 16:
        return 4;
      default:
        throw std::invalid_argument("tile_vectors must be 1, 2, 4, 8 or 16");
      }
    }
  }

  std::string active_isa_name(Isa isa) { return select(isa).name; }

  // Scheduling: points are grouped into tiles of `lanes = vectors * WIDTH` points, which one kernel call processes.
  // Tiles are either distributed over threads (default) or, when there are very few tiles, the threads work on
  // the m loop inside each tile (parallel_orders). The output array is indexed by the original point order.
  RawFields synth_points(const Prepared &prep, const Points &pts, int grad, const SynthOptions &opts)
  {
    if (grad < 0 || grad > 2)
      throw std::invalid_argument("synth_points: grad must be 0, 1 or 2");
    const std::int64_t np = pts.size();
    const int nf = kernels::fields_for_grad(grad);
    RawFields out;
    out.grad = grad;
    out.f.resize(nf);
    for (auto &v : out.f)
      v.assign(np, 0.0);
    if (np == 0)
      return out;

    const kernels::Table &tab = select(opts.isa);
    const int threads = opts.threads > 0 ? opts.threads : omp_get_max_threads();
    int vectors = opts.tile_vectors;
    // Default tile width (1, 2 or 4 vectors) chosen by measurement on AVX-512 and AVX2: the accumulators of the
    // kernels exhaust the vector registers at larger widths (measured in this project; the raw timings are not archived).
    if (vectors == 0)
      vectors = grad == 0 ? 4 : (grad == 1 ? 2 : 1);
    const int bi = b_index(vectors);
    const int lanes = vectors * tab.width;

    kernels::PrepView pv{prep.nmax, prep.rec->ab.data(), prep.cs.data(), prep.rec->f.data(), prep.rec->d.data(), prep.R};
    // Polar truncation threshold μ for the Jekeli-type rule: Jekeli, Lee and Kwon (2007), J. Geod. 81, 603-615.
    // The code uses μ = max(100, 0.01 * nmax). This choice is the project's; the benchmarks in benchmark/alf and
    // benchmark/meridian verify the resulting accuracy, no separate error bound was derived for it.
    const double pt = opts.polar_optimization ? std::max(100.0, 0.01 * prep.nmax) : -1.0;
    const kernels::TileOpts topts{pt, opts.skip_forbidden, opts.rotate_trig};

    // Sorting by |φ| makes the lanes of one tile share a similar forbidden (negligible) zone, so the skip test
    // in the kernel is effective. The output is still written back in the original order.
    std::vector<std::int64_t> order(np);
    std::iota(order.begin(), order.end(), std::int64_t{0});
    if (opts.sort_latitude)
      std::sort(order.begin(), order.end(), [&](std::int64_t a, std::int64_t b)
                { return std::fabs(pts.lat[a]) < std::fabs(pts.lat[b]); });

    const std::int64_t ntiles = (np + lanes - 1) / lanes;
    const bool by_order = opts.parallel_orders == 1 ||
                          (opts.parallel_orders < 0 && threads > 1 && ntiles < 2 * threads);

    const auto run_tile = [&](std::int64_t it, kernels::TileFn fn, double *lat, double *lon, double *r, double *buf)
    {
      const std::int64_t i0 = it * lanes;
      const int count = static_cast<int>(std::min<std::int64_t>(lanes, np - i0));
      for (int l = 0; l < lanes; ++l)
      {
        const std::int64_t id = order[i0 + std::min(l, count - 1)]; // padding: repeat the last point of the tile
        lat[l] = pts.lat[id];
        lon[l] = pts.lon[id];
        r[l] = pts.r[id];
      }
      fn(pv, topts, lat, lon, r, buf);
      for (int l = 0; l < count; ++l)
      {
        const std::int64_t id = order[i0 + l];
        for (int k = 0; k < nf; ++k)
          out.f[k][id] = buf[k * lanes + l];
      }
    };

    if (by_order)
    {
      // Few tiles: each tile is processed by all threads inside the m loop. Note that omp_set_num_threads changes the
      // process-wide OpenMP setting for later parallel regions (existing behaviour, kept unchanged).
      AlignedVec<double> lat(lanes), lon(lanes), r(lanes), buf(nf * lanes);
      omp_set_num_threads(threads);
      for (std::int64_t it = 0; it < ntiles; ++it)
        run_tile(it, tab.tile_mpar[grad][bi], lat.data(), lon.data(), r.data(), buf.data());
    }
    else
    {
      // Many tiles: one tile per work item, dynamic schedule because tiles can have different costs.
#pragma omp parallel num_threads(threads)
      {
        AlignedVec<double> lat(lanes), lon(lanes), r(lanes), buf(nf * lanes);
#pragma omp for schedule(dynamic, 1)
        for (std::int64_t it = 0; it < ntiles; ++it)
          run_tile(it, tab.tile[grad][bi], lat.data(), lon.data(), r.data(), buf.data());
      }
    }
    return out;
  }
}
