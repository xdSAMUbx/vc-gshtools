/* charm_grid_bench.c: CHarm driver of the grid synthesis benchmark (see benchmark/grid/code/README.md).
 *
 * Synthesises V and T on one latitude-longitude grid per run, with CHarm's grid-wise backward synthesis
 * (charm_shs_point() on a point structure of type GRID, GRID_GL, GRID_DH1 or GRID_DH2):
 *
 *   gl        Gauss-Legendre grid (Sneeuw, 1994): nlat = nmax + 1, nlon = 2 nmax + 2 (charm_crd_point_gl)
 *   dh1       non-equiangular Driscoll-Healy grid (Driscoll and Healy, 1994): 2(nmax+1) x 2(nmax+1)
 *             (charm_crd_point_dh1)
 *   dh2       equiangular Driscoll-Healy grid (Wieczorek and Meschede, 2018): 2(nmax+1) x 4(nmax+1)
 *             (charm_crd_point_dh2)
 *   regional  201 x 201 grid: latitude 45 -> 35 deg (north to south), longitude -110 -> -100 deg, step 0.05 deg,
 *             built with charm_crd_point_init(CHARM_CRD_POINT_GRID, ...)
 *
 * Grid latitudes are spherical (geocentric) latitudes, every point has the constant radius r = R of the model
 * header. Node coordinates are read back from the CHarm point structure, so the printed lat/lon are the ones
 * CHarm used.
 *
 * Definitions (benchmark/shs/CONTRACT.md, same as the point benchmark):
 *   V_m2s2 = charm_shs_point(): GM/r sum_n (R/r)^n sum_m (C cos m lon + S sin m lon) P_nm, r = R,
 *            degree 0 included (C00 as in the file). No normal potential is subtracted.
 *   T_m2s2 = V - U, with U the GRS80 normal potential at (r = R, geocentric latitude), closed form as in
 *            benchmark/shs/code/python/shs_common.py normal_potential_closed(). U depends on latitude only,
 *            so it is evaluated once per latitude row. The delta C00 term is not needed: the closed form U
 *            contains the GM term of the ellipsoid.
 *
 * Timing: the wall time (omp_get_wtime) covers charm_shs_point() only (V on the whole grid). U and T are
 * computed after the timed repetitions. Peak memory is the process peak RSS from getrusage(RUSAGE_SELF)
 * (ru_maxrss, kilobytes on Linux/WSL), so it is cumulative over the process and not per repetition.
 *
 * Output (never a whole grid, only the selected rows):
 *   <out-dir>/charm_<case>_rows.csv  header row,lat_deg,lon_deg,V_m2s2,T_m2s2 ; values with %.17g
 *   <out-dir>/timings.csv            appended, one line per repetition:
 *                                    engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb
 *   Rows: the index is 0-based and counts from the northernmost row. Global selection: the row whose latitude
 *   is closest to 0, 30, 60, 85 and 89.5 deg north and to -30, -60, -85 and -89.5 deg, plus rows 0 and nlat-1.
 *   Regional default: rows 0, 50, 100, 150 and 200 (45, 42.5, 40, 37.5 and 35 deg). --rows i,j,k overrides both.
 *
 * Parallelism: the OpenMP threads are set with omp_set_num_threads(--threads); CHarm's own loops (OpenMP) use them.
 * The FFT path of CHarm (FFTW) is selected by CHarm itself for the global grids; see recurrence_label().
 *
 * Usage:
 *   charm_grid_bench --gfc <tide_free.gfc> --nmax N --case gl|dh1|dh2|regional --threads K --repeat R
 *                    --out-dir <dir> [--rows i,j,k]
 *
 * Exit codes: 0 ok, 1 runtime or CHarm error, 2 bad arguments or model mismatch.
 *
 * Verified against the CHarm 0.4.11 sources on disk (build/external/charm): grid shapes and lon steps of
 * src/crd/crd_point_gl_chunk.c, crd_point_dh_chunk.c, crd_point_dh_lat_w_chunk.c and the *_shape.c files; the
 * charm_shs_point() contract in charm/charm_shs.h. NOT compiled or run by the author of this file (see README).
 *
 * References:
 *   Bucha, B. (2026) CHarm, Earth Science Informatics 19, 29.
 *   Sneeuw, N. (1994) Global spherical harmonic analysis by least-squares and numerical quadrature methods in
 *     historical perspective. Geophysical Journal International 118, 707-716.
 *   Driscoll, J. R., Healy, D. M. (1994) Computing Fourier transforms and convolutions on the 2-sphere.
 *     Advances in Applied Mathematics.
 *   Wieczorek, M. A., Meschede, M. (2018) SHTools: tools for working with spherical harmonics.
 *     Geochemistry, Geophysics, Geosystems.
 *   Balmino, G., Vales, N., Bonvalot, S., Briais, A. (2012) Spherical harmonic modelling to ultra-high degree of
 *     Bouguer and isostatic anomalies. Journal of Geodesy 86, 499-520.
 *   Heiskanen, W. A., Moritz, H. (1967) Physical Geodesy. W. H. Freeman, San Francisco.
 */
#include <charm/charm.h>

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>

#ifndef CHARM_VERSION_TAG
#  define CHARM_VERSION_TAG "0.4.11 (commit 1ead5828)"
#endif

/* GRS80 constants. Same literals as code/python/shs_common.py and benchmark/shs/code/c/charm_shs_bench.c. */
#define GRS80_A 6378137.0
#define GRS80_E2 0.00669438002290
#define GRS80_GM 3.986005e14
#define GRS80_OMEGA 7.292115e-5

#define PI_D 3.14159265358979323846
#define DEG2RAD (PI_D / 180.0)
#define RAD2DEG (180.0 / PI_D)

/* Regional grid (fixed by the benchmark definition). */
#define REG_LAT_NORTH 45.0
#define REG_LAT_SOUTH 35.0
#define REG_LON_WEST -110.0
#define REG_LON_EAST -100.0
#define REG_N 201 /* points per direction, step 0.05 deg inclusive */

#define MAX_SEL 64

#define CHECK(err)                                   \
    do                                               \
    {                                                \
        if (!charm_err_isempty(err))                 \
            charm_err_handler(err, 1);               \
    } while (0)

typedef struct
{
    const char *gfc;
    long nmax;
    const char *grid_case;
    int threads;
    int repeat;
    const char *out_dir;
    const char *rows_arg; /* NULL when not given */
} Args;

static void die(const char *msg)
{
    fprintf(stderr, "charm_grid_bench: %s\n", msg);
    exit(2);
}

static void usage(void)
{
    fprintf(stderr,
            "usage: charm_grid_bench --gfc F --nmax N --case gl|dh1|dh2|regional --threads K --repeat R\n"
            "                        --out-dir DIR [--rows i,j,k]\n");
    exit(2);
}

static Args parse_args(int argc, char **argv)
{
    Args a;
    int have_gfc = 0, have_nmax = 0, have_case = 0, have_out = 0;
    a.gfc = NULL;
    a.nmax = -1;
    a.grid_case = NULL;
    a.threads = 1;
    a.repeat = 1;
    a.out_dir = NULL;
    a.rows_arg = NULL;
    for (int k = 1; k < argc; k++)
    {
        const char *arg = argv[k];
        if (k + 1 >= argc && strcmp(arg, "--help") != 0)
        {
            fprintf(stderr, "charm_grid_bench: missing value for %s\n", arg);
            usage();
        }
        if (strcmp(arg, "--gfc") == 0) { a.gfc = argv[++k]; have_gfc = 1; }
        else if (strcmp(arg, "--nmax") == 0) { a.nmax = strtol(argv[++k], NULL, 10); have_nmax = 1; }
        else if (strcmp(arg, "--case") == 0) { a.grid_case = argv[++k]; have_case = 1; }
        else if (strcmp(arg, "--threads") == 0) { a.threads = (int)strtol(argv[++k], NULL, 10); }
        else if (strcmp(arg, "--repeat") == 0) { a.repeat = (int)strtol(argv[++k], NULL, 10); }
        else if (strcmp(arg, "--out-dir") == 0) { a.out_dir = argv[++k]; have_out = 1; }
        else if (strcmp(arg, "--rows") == 0) { a.rows_arg = argv[++k]; }
        else if (strcmp(arg, "--help") == 0) usage();
        else
        {
            fprintf(stderr, "charm_grid_bench: unknown argument %s\n", arg);
            usage();
        }
    }
    if (!(have_gfc && have_nmax && have_case && have_out))
        usage();
    if (a.nmax < 0 || a.threads < 1 || a.repeat < 1)
        die("nmax must be >= 0, threads and repeat >= 1");
    if (strcmp(a.grid_case, "gl") != 0 && strcmp(a.grid_case, "dh1") != 0 && strcmp(a.grid_case, "dh2") != 0 &&
        strcmp(a.grid_case, "regional") != 0)
        die("--case must be gl, dh1, dh2 or regional");
    return a;
}

/* Wall clock is omp_get_wtime(); peak memory is the process peak RSS (ru_maxrss, KB on Linux). */
static double peak_mem_mb(void)
{
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0)
        return -1.0;
    return (double)ru.ru_maxrss / 1024.0;
}

/* Header keys used: earth_gravity_constant (or gravity_constant), radius, max_degree. Same keys as the SHS driver. */
static void read_header(const char *path, double *gm, double *radius, long *max_degree)
{
    FILE *f = fopen(path, "r");
    char line[1024], key[128], val[256];
    int have_gm = 0, have_r = 0, have_max = 0, ended = 0;
    if (!f)
        die("cannot open the .gfc file");
    while (fgets(line, sizeof line, f))
    {
        if (strncmp(line, "end_of_head", 11) == 0)
        {
            ended = 1;
            break;
        }
        if (sscanf(line, "%127s %255s", key, val) != 2)
            continue;
        if (strcmp(key, "earth_gravity_constant") == 0 || strcmp(key, "gravity_constant") == 0)
        {
            *gm = strtod(val, NULL);
            have_gm = 1;
        }
        else if (strcmp(key, "radius") == 0)
        {
            *radius = strtod(val, NULL);
            have_r = 1;
        }
        else if (strcmp(key, "max_degree") == 0)
        {
            *max_degree = strtol(val, NULL, 10);
            have_max = 1;
        }
    }
    fclose(f);
    if (!ended)
        die("no end_of_head line in the model file");
    if (!(have_gm && have_r && have_max))
        die("model header lacks earth_gravity_constant, radius or max_degree");
}

/* q(x) of normal_potential_closed() in shs_common.py */
static double q_func(double x, double e_lin)
{
    return 0.5 * ((1.0 + 3.0 * x * x / (e_lin * e_lin)) * atan(e_lin / x) - 3.0 * x / e_lin);
}

/* normal_potential_closed() of shs_common.py: gravitational part of the GRS80 normal potential [m^2/s^2] */
static double normal_potential_closed(double r, double phi_c)
{
    double e_lin = GRS80_A * sqrt(GRS80_E2);
    double b_semi = GRS80_A * sqrt(1.0 - GRS80_E2);
    double z = r * sin(phi_c);
    double d = r * r - e_lin * e_lin;
    double u2 = 0.5 * (d + sqrt(d * d + 4.0 * e_lin * e_lin * z * z));
    double u = sqrt(u2);
    double sin2_beta = z * z / u2;
    double q0 = q_func(b_semi, e_lin);
    double flattening = 0.5 * GRS80_OMEGA * GRS80_OMEGA * GRS80_A * GRS80_A * (q_func(u, e_lin) / q0) *
                        (sin2_beta - 1.0 / 3.0);
    return GRS80_GM / e_lin * atan(e_lin / u) + flattening;
}

/* Recurrence label as CHarm documents its choice (charm_shs.h, charm_shs_point): the FFT path is used on a grid
 * when nlon > 1, (nlon - 1) / 2 >= nmax (integer division), constant dlon, lon[0] == 0 and
 * lon[nlon-1] + dlon == 2 pi; otherwise Chebyshev recurrences. The label is derived here from those documented
 * conditions. It is NOT read back from CHarm, so it is an expectation, not a measurement. */
static const char *recurrence_label(const charm_point *pnt, long nmax)
{
    const double thr = 1e-9;
    if (pnt->nlon < 2)
        return "chebyshev";
    double dlon = pnt->lon[1] - pnt->lon[0];
    long nlon = (long)pnt->nlon;
    int ok = ((nlon - 1) / 2 >= nmax) && (fabs(pnt->lon[0]) < thr) &&
             (fabs(pnt->lon[nlon - 1] + dlon - 2.0 * PI_D) < thr);
    for (long j = 2; j < nlon && ok; j++)
        if (fabs((pnt->lon[j] - pnt->lon[j - 1]) - dlon) > thr)
            ok = 0;
    return ok ? "fft" : "chebyshev";
}

static void append_timing(const char *out_dir, const char *grid_case, const char *recurrence, int threads,
                          long nmax, size_t nlat, size_t nlon, int rep, double wall, double peak)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/timings.csv", out_dir);
    FILE *probe = fopen(path, "r");
    int need_header = (probe == NULL);
    if (probe)
        fclose(probe);
    FILE *f = fopen(path, "a");
    if (!f)
        die("cannot append to timings.csv (does --out-dir exist?)");
    if (need_header)
        fprintf(f, "engine,case,recurrence,threads,nmax,nlat,nlon,rep,wall_s,peak_mem_mb\n");
    fprintf(f, "charm,%s,%s,%d,%ld,%zu,%zu,%d,%.9g,%.3f\n", grid_case, recurrence, threads, nmax, nlat, nlon, rep,
            wall, peak);
    fclose(f);
}

/* Index of the row whose latitude is closest to target_deg (lat in radians, row 0 = northernmost). */
static long nearest_row(const double *lat_rad, long nlat, double target_deg)
{
    long best = 0;
    double bd = 1e300;
    for (long i = 0; i < nlat; i++)
    {
        double d = fabs(lat_rad[i] * RAD2DEG - target_deg);
        if (d < bd)
        {
            bd = d;
            best = i;
        }
    }
    return best;
}

static void add_row(long *sel, int *nsel, long v)
{
    for (int k = 0; k < *nsel; k++)
        if (sel[k] == v)
            return;
    if (*nsel >= MAX_SEL)
        die("too many selected rows");
    sel[(*nsel)++] = v;
}

static void sort_rows(long *sel, int n)
{
    for (int i = 1; i < n; i++)
        for (int k = i; k > 0 && sel[k - 1] > sel[k]; k--)
        {
            long t = sel[k - 1];
            sel[k - 1] = sel[k];
            sel[k] = t;
        }
}

/* Selected rows: the default of the case, or the list of --rows (validated against nlat). */
static int select_rows(const Args *a, const double *lat_rad, long nlat, long *sel)
{
    int nsel = 0;
    if (a->rows_arg)
    {
        const char *p = a->rows_arg;
        while (*p)
        {
            char *end = NULL;
            long v = strtol(p, &end, 10);
            if (end == p || v < 0 || v >= nlat)
                die("--rows must be a comma-separated list of row indices in [0, nlat-1]");
            add_row(sel, &nsel, v);
            p = end;
            if (*p == ',')
                p++;
            else if (*p != '\0')
                die("--rows must be a comma-separated list of row indices");
        }
    }
    else if (strcmp(a->grid_case, "regional") == 0)
    {
        add_row(sel, &nsel, 0);
        add_row(sel, &nsel, (nlat - 1) / 4);
        add_row(sel, &nsel, (nlat - 1) / 2);
        add_row(sel, &nsel, 3 * (nlat - 1) / 4);
        add_row(sel, &nsel, nlat - 1);
    }
    else
    {
        static const double targets[] = {0.0, 30.0, 60.0, 85.0, 89.5, -30.0, -60.0, -85.0, -89.5};
        for (int k = 0; k < (int)(sizeof targets / sizeof targets[0]); k++)
            add_row(sel, &nsel, nearest_row(lat_rad, nlat, targets[k]));
        add_row(sel, &nsel, 0);
        add_row(sel, &nsel, nlat - 1);
    }
    sort_rows(sel, nsel);
    return nsel;
}

static void write_rows(const Args *a, const charm_point *pnt, const long *sel, int nsel, const double *f,
                       const double *u_row, long nlat, long nlon)
{
    char path[4096];
    snprintf(path, sizeof path, "%s/charm_%s_rows.csv", a->out_dir, a->grid_case);
    FILE *out = fopen(path, "w");
    if (!out)
        die("cannot write the rows CSV (does --out-dir exist?)");
    fprintf(out, "row,lat_deg,lon_deg,V_m2s2,T_m2s2\n");
    for (int k = 0; k < nsel; k++)
    {
        long i = sel[k];
        double lat_deg = pnt->lat[i] * RAD2DEG;
        for (long j = 0; j < nlon; j++)
        {
            double v = f[i * nlon + j];
            double t = v - u_row[i];
            fprintf(out, "%ld,%.17g,%.17g,%.17g,%.17g\n", i, lat_deg, pnt->lon[j] * RAD2DEG, v, t);
        }
    }
    fclose(out);
    (void)nlat;
}

int main(int argc, char **argv)
{
    Args a = parse_args(argc, argv);
    int regional = (strcmp(a.grid_case, "regional") == 0);

    double gm = 0.0, radius = 0.0;
    long max_degree = 0;
    read_header(a.gfc, &gm, &radius, &max_degree);
    if (a.nmax > max_degree)
    {
        fprintf(stderr, "charm_grid_bench: nmax=%ld exceeds max_degree=%ld of the model\n", a.nmax, max_degree);
        return 2;
    }

    /* Nested OpenMP off: CHarm's loops use the threads set here, the driver adds no parallel region. */
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);
    omp_set_num_threads(a.threads);

    charm_err *err = charm_err_init();
    if (!err)
        die("charm_err_init failed");
    charm_shc *shcs = charm_shc_calloc((unsigned long)a.nmax, 1.0, 1.0);
    if (!shcs)
        die("charm_shc_calloc failed");
    charm_shc_read_gfc(a.gfc, (unsigned long)a.nmax, NULL, shcs, err);
    CHECK(err);
    if (fabs(shcs->mu - gm) > 1e-9 * gm || fabs(shcs->r - radius) > 1e-6)
    {
        fprintf(stderr, "charm_grid_bench: GM/R mismatch: CHarm (%.10e, %.6f) vs header (%.10e, %.6f)\n", shcs->mu,
                shcs->r, gm, radius);
        return 2;
    }

    /* ---- grid ---- */
    charm_point *pnt = NULL;
    double *reg_lat = NULL, *reg_lon = NULL, *reg_r = NULL;
    if (strcmp(a.grid_case, "gl") == 0)
        pnt = charm_crd_point_gl((unsigned long)a.nmax, radius);
    else if (strcmp(a.grid_case, "dh1") == 0)
        pnt = charm_crd_point_dh1((unsigned long)a.nmax, radius);
    else if (strcmp(a.grid_case, "dh2") == 0)
        pnt = charm_crd_point_dh2((unsigned long)a.nmax, radius);
    else
    {
        reg_lat = malloc(REG_N * sizeof *reg_lat);
        reg_lon = malloc(REG_N * sizeof *reg_lon);
        reg_r = malloc(REG_N * sizeof *reg_r);
        if (!reg_lat || !reg_lon || !reg_r)
            die("allocation failed (regional grid)");
        for (int i = 0; i < REG_N; i++)
        {
            /* 45 -> 35 deg, north to south, exact end points */
            double lat_deg = REG_LAT_NORTH - (REG_LAT_NORTH - REG_LAT_SOUTH) * i / (double)(REG_N - 1);
            reg_lat[i] = lat_deg * DEG2RAD; /* spherical (geocentric) latitude */
            reg_r[i] = radius;
        }
        for (int j = 0; j < REG_N; j++)
            reg_lon[j] = (REG_LON_WEST + (REG_LON_EAST - REG_LON_WEST) * j / (double)(REG_N - 1)) * DEG2RAD;
        pnt = charm_crd_point_init(CHARM_CRD_POINT_GRID, REG_N, REG_N, reg_lat, reg_lon, reg_r);
    }
    if (!pnt)
        die("CHarm grid construction failed");

    long nlat = (long)pnt->nlat, nlon = (long)pnt->nlon;
    if (!(pnt->lat[0] > pnt->lat[nlat - 1]))
        die("CHarm grid is not ordered north to south (row 0 must be the northernmost)");
    const char *recurrence = recurrence_label(pnt, a.nmax);

    double *f = malloc((size_t)pnt->npoint * sizeof *f);
    double *u_row = malloc((size_t)nlat * sizeof *u_row);
    if (!f || !u_row)
        die("allocation failed (values)");

    /* ---- timed repetitions: V on the whole grid ---- */
    for (int rep = 1; rep <= a.repeat; rep++)
    {
        double t0 = omp_get_wtime();
        charm_shs_point(pnt, shcs, (unsigned long)a.nmax, f, err);
        double wall = omp_get_wtime() - t0;
        CHECK(err);
        append_timing(a.out_dir, a.grid_case, recurrence, a.threads, a.nmax, (size_t)nlat, (size_t)nlon, rep, wall,
                      peak_mem_mb());
        fprintf(stdout, "repeat %d: wall=%.6f s, %s\n", rep, wall, recurrence);
    }

    /* ---- T = V - U per latitude row, then the selected rows ---- */
    for (long i = 0; i < nlat; i++)
        u_row[i] = normal_potential_closed(radius, pnt->lat[i]);

    long *sel = malloc(MAX_SEL * sizeof *sel);
    if (!sel)
        die("allocation failed (rows)");
    int nsel = select_rows(&a, pnt->lat, nlat, sel);
    write_rows(&a, pnt, sel, nsel, f, u_row, nlat, nlon);

    fprintf(stdout, "charm_grid_bench: %s nmax=%ld nlat=%ld nlon=%ld recurrence=%s CHarm %s, rows written: %d\n",
            a.grid_case, a.nmax, nlat, nlon, recurrence, CHARM_VERSION_TAG, nsel);

    free(sel);
    free(u_row);
    free(f);
    if (regional)
    {
        free(reg_lat);
        free(reg_lon);
        free(reg_r);
    }
    charm_crd_point_free(pnt);
    charm_shc_free(shcs);
    charm_err_free(err);
    return 0;
}
