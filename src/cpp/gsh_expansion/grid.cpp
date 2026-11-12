#include "grid.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <omp.h>
#include <stdexcept>
#include "kernels.hpp"

namespace gsh
{
  using namespace kernels;
  namespace
  {
    using Clock = std::chrono::steady_clock;
    double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

    // Kernel table of the requested ISA (Auto selects the best available one, see shs.cpp).
    const Table &auto_table(Isa isa)
    {
      const std::string n = active_isa_name(isa);
      return n == "avx512" ? table_avx512() : (n == "avx2" ? table_avx2() : table_scalar());
    }

    // Gauss-Legendre nodes: Newton iteration on P_L from the standard initial guess cos(π(i + 3/4)/(L + 1/2)),
    // then exact symmetrisation so that the nodes come in exact ± pairs.
    void gauss_nodes(int L, AlignedVec<double> &lat)
    {
      lat.resize(L);
#pragma omp parallel for schedule(static)
      for (int i = 0; i < L; ++i)
      {
        double z = std::cos(PI * (i + 0.75) / (L + 0.5));
        for (int it = 0; it < 100; ++it)
        {
          double p0 = 1, p1 = z;
          for (int k = 2; k <= L; ++k)
          {
            const double p2 = ((2.0 * k - 1) * z * p1 - (k - 1.0) * p0) / k;
            p0 = p1;
            p1 = p2;
          }
          const double dp = L * (z * p1 - p0) / (z * z - 1);
          const double dz = p1 / dp;
          z -= dz;
          if (std::abs(dz) < 1e-16)
            break;
        }
        lat[i] = std::asin(z);
      }
      // Exact antisymmetry: lat[L-1-i] = -lat[i]; the middle node of an odd-length set is exactly 0.
      for (int i = 0; i < L / 2; ++i)
      {
        const double v = 0.5 * (lat[i] - lat[L - 1 - i]);
        lat[i] = v;
        lat[L - 1 - i] = -v;
      }
      if (L % 2)
        lat[L / 2] = 0.0;
    }

    // Global grid with one radius for all rows and uniform longitudes starting at 0.
    Grid uniform(AlignedVec<double> lat, std::int64_t nlon, double r)
    {
      Grid g;
      g.r.assign(lat.size(), r);
      g.lat = std::move(lat);
      g.nlon = nlon;
      g.dlon = 2.0 * PI / static_cast<double>(nlon);
      return g;
    }

    // Which lump group feeds each field, and the m-weight it carries: kind 0 = plain (a, b); kind 1 = m (b, -a),
    // which is the longitude derivative; kind 2 = m² (a, b), the second longitude derivative. See synth_grid.
    struct FieldMap
    {
      int group, kind;
    };
    constexpr FieldMap field_map[F_COUNT] = {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {3, 0}, {4, 0}, {1, 1}, {5, 0}, {2, 1}, {0, 2}};
    // Parity under φ -> -φ. P̄_nm(-φ) = (-1)^(n-m) P̄_nm(φ), so terms with even n - m keep their sign and odd ones flip.
    // Groups that contain a φ-derivative (dP) pick up one extra sign flip. group_is_dp marks those groups.
    constexpr bool group_is_dp[6] = {false, false, true, false, true, false};

    struct LatTask
    {
      double lat, r;
      std::int64_t row, row_mirror; // row_mirror = -1 when there is no opposite latitude
    };

    // Sorts rows by |φ| and pairs each latitude +φ with -φ when both exist with the same radius. A mirror pair is
    // computed once; its partner row is produced with the parity signs (see emit). Rows with φ = 0 are never paired.
    std::vector<LatTask> build_tasks(const Grid &g)
    {
      const std::int64_t n = g.nlat();
      std::vector<std::int64_t> id(n);
      for (std::int64_t i = 0; i < n; ++i)
        id[i] = i;
      std::sort(id.begin(), id.end(), [&](auto a, auto b)
                {
        const double x = std::abs(g.lat[a]), y = std::abs(g.lat[b]);
        return x != y ? x < y : g.lat[a] > g.lat[b]; });
      std::vector<LatTask> t;
      for (std::int64_t k = 0; k < n; ++k)
      {
        const auto i = id[k];
        LatTask task{g.lat[i], g.r[i], i, -1};
        if (g.lat[i] != 0.0 && k + 1 < n)
        {
          const auto j = id[k + 1];
          if (g.lat[j] == -g.lat[i] && g.r[j] == g.r[i])
          {
            task.row_mirror = j;
            ++k;
          }
        }
        t.push_back(task);
      }
      return t;
    }

    // Grid lumps from the full MSFC triangle, for rows near the poles where the column-wise derivative (1/cos φ form)
    // is inaccurate. The output layout is the one of grid_lumps (column_kernel.inl), for a single lane.
    void oracle_lumps(const Prepared &prep, int grad, double lat, double r, double *L, int lanes, int lane)
    {
      const int nmax = prep.nmax, P = grid_groups(grad);
      static_cast<void>(lanes);
      const bool pole = std::abs(lat) == PI / 2;
      const double t = pole ? (lat > 0 ? 1.0 : -1.0) : std::sin(lat), u = pole ? 0.0 : std::cos(lat);
      legendre::Coefs c;
      c.nmax = nmax;
      c.a.resize(tri_size(nmax));
      c.b.resize(tri_size(nmax));
      c.d = prep.rec->d;
      for (std::int64_t i = 0; i < tri_size(nmax); ++i)
      {
        c.a[i] = prep.rec->ab[2 * i];
        c.b[i] = prep.rec->ab[2 * i + 1];
      }
      std::vector<double> Pn(tri_size(nmax)), dP(tri_size(nmax)), ddP(tri_size(nmax));
      legendre::msfc(t, u, c, Pn);
      legendre::dphi(nmax, Pn, dP);
      legendre::dphi(nmax, dP, ddP);
      const double q = prep.R / r;
      for (int m = 0; m <= nmax; ++m)
      {
        const auto off = col_offset(m, nmax);
        double qn = std::pow(q, m);
        for (int n = m; n <= nmax; ++n, qn *= q)
        {
          const auto i = off + (n - m);
          const double C = prep.cs[2 * i], S = prep.cs[2 * i + 1], N = n;
          const int p = (n - m) & 1;
          const std::int64_t base = (static_cast<std::int64_t>(m) * 2 + p) * 2 * P;
          // Group values, in the same order as grid_lumps: P, (n+1)P, dP, (n+1)(n+2)P, (n+1)dP, ddP, all times q^n.
          const double vals[6] = {qn * Pn[i], qn * (N + 1) * Pn[i], qn * dP[i], qn * (N + 1) * (N + 2) * Pn[i], qn * (N + 1) * dP[i], qn * ddP[i]};
          for (int g = 0; g < P; ++g)
          {
            L[(base + 2 * g) * lanes + lane] += vals[g] * C;
            L[(base + 2 * g + 1) * lanes + lane] += vals[g] * S;
          }
        }
      }
    }
  }

  Grid Grid::gauss_legendre(int nmax, double r)
  {
    const int L = nmax + 1;
    AlignedVec<double> lat;
    gauss_nodes(L, lat);
    return uniform(std::move(lat), 2 * static_cast<std::int64_t>(L), r);
  }

  Grid Grid::driscoll_healy1(int nmax, double r)
  {
    const int L = nmax + 1;
    AlignedVec<double> lat(2 * L);
    for (int i = 0; i < 2 * L; ++i)
      lat[i] = i == L ? 0.0 : PI / 2 - i * PI / (2 * L);
    return uniform(std::move(lat), 2 * static_cast<std::int64_t>(L), r);
  }

  Grid Grid::driscoll_healy2(int nmax, double r)
  {
    Grid g = driscoll_healy1(nmax, r);
    g.nlon *= 2;
    g.dlon /= 2;
    return g;
  }

  Grid Grid::regular(double lat_north, double lat_south, std::int64_t nlat, double lon_west, double lon_east, std::int64_t nlon, double r)
  {
    Grid g;
    g.lat.resize(nlat);
    for (std::int64_t i = 0; i < nlat; ++i)
      g.lat[i] = deg2rad(nlat > 1 ? lat_north + (lat_south - lat_north) * i / (nlat - 1) : lat_north);
    g.r.assign(nlat, r);
    g.nlon = nlon;
    g.lon0 = deg2rad(lon_west);
    g.dlon = nlon > 1 ? deg2rad((lon_east - lon_west) / (nlon - 1)) : 0.0;
    return g;
  }

  bool Grid::full_circle() const
  {
    return nlon >= 2 && std::abs(dlon * nlon - 2.0 * PI) < 1e-12;
  }

  Points grid_to_points(const Grid &g)
  {
    Points p;
    const auto total = g.nlat() * g.nlon;
    p.lat.resize(total);
    p.lon.resize(total);
    p.r.resize(total);
    for (std::int64_t i = 0; i < g.nlat(); ++i)
      for (std::int64_t j = 0; j < g.nlon; ++j)
      {
        const auto k = i * g.nlon + j;
        p.lat[k] = g.lat[i];
        p.lon[k] = g.lon0 + j * g.dlon;
        p.r[k] = g.r[i];
      }
    return p;
  }

  RawFields synth_grid(const Prepared &prep, const Grid &grid, int grad, const SynthOptions &opts, GridTiming *timing)
  {
    if (grad < 0 || grad > 2)
      throw std::invalid_argument("synth_grid: grad must be 0, 1 or 2");
    const auto t_all = Clock::now();
    const int nf = fields_for_grad(grad), P = grid_groups(grad), nmax = prep.nmax;
    const Table &tab = auto_table(opts.isa);
    const int lanes = tab.width;
    const std::int64_t nlat = grid.nlat(), nlon = grid.nlon;
    if (grid.r.size() != grid.lat.size())
      throw std::invalid_argument("synth_grid: one radius per latitude is required");

    RawFields out;
    out.grad = grad;
    out.f.resize(nf);
    for (auto &v : out.f)
      v.resize(nlat * nlon);

    const bool use_fft = grid.full_circle();
    const auto tasks_all = build_tasks(grid);

    // With derivatives, rows with cos φ < 5e-3 are sent to the full-triangle oracle (same threshold as
    // FunctionalRequest::pole_cos_threshold). The value 5e-3 is a choice of this project.
    std::vector<LatTask> tasks, polar;
    for (const auto &t : tasks_all)
      (grad > 0 && std::cos(t.lat) < 5e-3 ? polar : tasks).push_back(t);

    PrepView pv{nmax, prep.rec->ab.data(), prep.cs.data(), prep.rec->f.data(), prep.rec->d.data(), prep.R};
    const double pt = opts.polar_optimization ? std::max(100.0, 0.01 * nmax) : -1.0;
    const TileOpts topts{pt, opts.skip_forbidden, opts.rotate_trig};

    // One complex-to-real inverse FFT of length nlon per field and row. FFTW_ESTIMATE avoids planning measurements;
    // FFTW_DESTROY_INPUT is allowed because the spectrum buffer is rebuilt for every transform.
    fftw_plan plan = nullptr;
    if (use_fft)
    {
      AlignedVec<double> tmp(nlon + 2);
      plan = fftw_plan_dft_c2r_1d(static_cast<int>(nlon), reinterpret_cast<fftw_complex *>(tmp.data()), tmp.data(), FFTW_ESTIMATE | FFTW_DESTROY_INPUT | FFTW_UNALIGNED);
      if (!plan)
        throw std::runtime_error("FFTW: could not create the plan");
    }

    // Turns the lumps of one row into one output row per field. For each field and order m it builds the
    // coefficients (a, b) of cos(mλ) and sin(mλ) from the even and odd lump pairs, then sums over m either
    // in the Fourier domain (FFT path) or directly (regional path). mirror = true emits the row at -φ.
    const auto emit = [&](const double *lump, int lane, std::int64_t row, bool mirror, std::vector<std::complex<double>> &spec, AlignedVec<double> &cj, AlignedVec<double> &sj)
    {
      const std::int64_t per_m = 2 * 2 * P * lanes;
      for (int field = 0; field < nf; ++field)
      {
        const FieldMap map = field_map[field];
        const bool dp = group_is_dp[map.group];
        // Signs for the even (se) and odd (so) (n - m) terms when the row is mirrored to -φ. Derivative groups get an
        // extra flip (see group_is_dp).
        const double se = mirror ? (dp ? -1.0 : 1.0) : 1.0;
        const double so = mirror ? (dp ? 1.0 : -1.0) : 1.0;
        double *const dst = out.f[field].data() + row * nlon;
        const auto ab_of = [&](int m, double &a, double &b)
        {
          const double *L = lump + m * per_m;
          const int l = 2 * map.group;
          const double ae = L[l * lanes + lane], be = L[(l + 1) * lanes + lane];
          const double ao = L[(2 * P + l) * lanes + lane], bo = L[(2 * P + l + 1) * lanes + lane];
          const double A = se * ae + so * ao, Bv = se * be + so * bo;
          const double M = m;
          if (map.kind == 0)
          {
            a = A;
            b = Bv;
          }
          else if (map.kind == 1)
          {
            // m (b cos - a sin) written as a' cos + b' sin with a' = m b and b' = -m a.
            a = M * Bv;
            b = -M * A;
          }
          else
          {
            a = M * M * A;
            b = M * M * Bv;
          }
        };

        if (use_fft)
        {
          // FFT path. The sum over m of Re[(a - i b) e^{i m λ0} e^{i m j Δ}] is written as a half-complex spectrum.
          // Frequency m is folded onto k = m mod nlon in [0, nlon/2] (aliasing, which is exact because the grid is
          // uniform). Each k in (0, nlon/2) gets half of the complex term, so that the inverse real transform
          // reproduces Re of the sum; k = 0 and the Nyquist bin keep only the real part.
          std::fill(spec.begin(), spec.end(), std::complex<double>(0.0, 0.0));
          for (int m = 0; m <= nmax; ++m)
          {
            double a, b;
            ab_of(m, a, b);
            if (a == 0.0 && b == 0.0)
              continue;
            // Re[(a - i b) e^{i m λ0} e^{i m j Δ}]
            std::complex<double> z(a, -b);
            if (grid.lon0 != 0.0)
              z *= std::polar(1.0, m * grid.lon0);
            std::int64_t k = m % nlon;
            if (k > nlon - k)
            {
              k = nlon - k;
              z = std::conj(z);
            }
            if (k == 0 || 2 * k == nlon)
              spec[k] += std::complex<double>(z.real(), 0.0);
            else
              spec[k] += 0.5 * z;
          }
          fftw_execute_dft_c2r(plan, reinterpret_cast<fftw_complex *>(spec.data()), dst);
        }
        else
        {
          // Direct path (regional or non-uniform grids). The phases e^{i m λ_j} are advanced by the rotation
          // (c, s) -> (c c1 - s s1, s c1 + c s1) with the step angle Δλ = dlon.
          std::fill(dst, dst + nlon, 0.0);
          for (int m = 0; m <= nmax; ++m)
          {
            double a, b;
            ab_of(m, a, b);
            if (a == 0.0 && b == 0.0)
              continue;
            const double c1 = std::cos(m * grid.dlon), s1 = std::sin(m * grid.dlon);
            double c = std::cos(m * grid.lon0), s = std::sin(m * grid.lon0);
            for (std::int64_t j = 0; j < nlon; ++j)
            {
              dst[j] += a * c + b * s;
              const double c2 = c * c1 - s * s1;
              s = s * c1 + c * s1;
              c = c2;
            }
          }
        }
      }
    };

    const std::int64_t ntiles = (static_cast<std::int64_t>(tasks.size()) + lanes - 1) / lanes;
    const int threads = opts.threads > 0 ? opts.threads : omp_get_max_threads();
    double t_lump = 0, t_fft = 0;
    const std::size_t lump_size = static_cast<std::size_t>(nmax + 1) * 2 * 2 * P * lanes;

    // Tiles of `lanes` rows are distributed over threads. Each thread owns its lump buffer and its spectrum, so
    // there is no sharing; the timers are reduced across threads (thread-seconds).
#pragma omp parallel num_threads(threads) reduction(+ : t_lump, t_fft)
    {
      AlignedVec<double> lump(lump_size), lat(lanes), r(lanes), cj, sj;
      std::vector<std::complex<double>> spec(nlon / 2 + 1);
#pragma omp for schedule(dynamic, 1)
      for (std::int64_t it = 0; it < ntiles; ++it)
      {
        const std::int64_t i0 = it * lanes;
        const int count = static_cast<int>(std::min<std::int64_t>(lanes, static_cast<std::int64_t>(tasks.size()) - i0));
        for (int l = 0; l < lanes; ++l)
        {
          const auto &t = tasks[i0 + std::min(l, count - 1)];
          lat[l] = t.lat;
          r[l] = t.r;
        }
        auto t0 = Clock::now();
        tab.grid[grad](pv, topts, 0, nmax, lat.data(), r.data(), lump.data());
        t_lump += secs(t0);
        t0 = Clock::now();
        for (int l = 0; l < count; ++l)
        {
          const auto &t = tasks[i0 + l];
          emit(lump.data(), l, t.row, false, spec, cj, sj);
          if (t.row_mirror >= 0)
            emit(lump.data(), l, t.row_mirror, true, spec, cj, sj);
        }
        t_fft += secs(t0);
      }
    }

    // Near-polar rows: full triangle per row (slow, but there are few such rows). This loop is serial.
    for (const auto &t : polar)
    {
      AlignedVec<double> lump(lump_size, 0.0), cj, sj;
      std::vector<std::complex<double>> spec(nlon / 2 + 1);
      oracle_lumps(prep, grad, t.lat, t.r, lump.data(), lanes, 0);
      emit(lump.data(), 0, t.row, false, spec, cj, sj);
      if (t.row_mirror >= 0)
        emit(lump.data(), 0, t.row_mirror, true, spec, cj, sj);
    }

    if (plan)
      fftw_destroy_plan(plan);
    if (timing)
    {
      timing->lumps = t_lump;
      timing->fft = t_fft;
      timing->total = secs(t_all);
    }
    return out;
  }
}
