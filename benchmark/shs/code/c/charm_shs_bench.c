/* charm_shs_bench.c: driver of the CHarm engine for the spherical-harmonic synthesis benchmark.
 *
 * Definitions: benchmark/shs/results/csv.md. Identical to code/python/run_pyharm.py and
 * code/python/shs_common.py (the python functions are ported line by line below, same formulas, same order of
 * operations where practical). One run evaluates at each point, in ONE synthesis pass, the three functionals:
 *   V_m2s2   = charm_shs_point()                       gravitational potential of the model, mu/r sum,
 *                                                       degree 0 included (C00 as in the file)
 *   T_m2s2   = V - U                                   U: GRS80 level-ellipsoid potential, closed form
 *   dg_mgal  = (|g| - gamma_SH) * 1e5                  g from charm_shs_point_grad1() (LNOF), with the
 *                                                       centrifugal term; gamma_SH: project series
 * Geodetic (GRS80) -> geocentric (phi', lambda, r) conversion is exact and identical to the python runner.
 *
 * Parallelisation (see benchmark/shs/code/c/README.md), selected by --parallel:
 *   * charm (default, the official form used by run_all.ps1): one CHarm call over all points (as pyharm calls
 *     CHarm); CHarm's OpenMP over m uses K threads and the driver adds no parallelism.
 *   * points (DIAGNOSTIC, not in run_all.ps1): the driver parallelises over POINTS with OpenMP (--threads K). Each point is a one-point
 *     CHarm point structure (charm_crd_point_init, borrowing the driver's arrays) with its own call. CHarm's
 *     own scattered-point synthesis is OpenMP-parallel over the spherical-harmonic orders m; inside the driver
 *     this inner region is kept at one thread, so the total number of threads never exceeds K.
 *   CHarm's FFT (FFTW) is not used on either path (it is for grids only).
 *
 * Usage:
 *   charm_shs_bench --gfc <tide_free.gfc> --points <csv> --nmax N --out <result.csv>
 *                   --timings <timings.csv> --set colorado|latsweep --threads K --repeat R
 *                   [--parallel charm|points]   (default charm)
 * Writes <out> (id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal; 15 significant digits), the companion <out>.md
 * and appends one row per repetition to <timings> (header written when the file is new or empty).
 * Wall time covers the synthesis only (all points, V and gradient); model loading is reported as load_s in notes.
 *
 * All text in English. */

#include <charm/charm.h>

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#endif

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/* Build flags as recorded in the notes; the build script overrides this macro. */
#ifndef CHARM_DRIVER_FLAGS
#  define CHARM_DRIVER_FLAGS "unknown"
#endif
/* CHarm version and commit (tag v0.4.11). The build script overrides this macro when needed. */
#ifndef CHARM_VERSION_TAG
#  define CHARM_VERSION_TAG "0.4.11 (commit 1ead5828)"
#endif

/* GRS80 constants. Same literals as code/python/shs_common.py (and src/cpp/gsh_utils/ellipsoid.hpp). */
#define GRS80_A 6378137.0
#define GRS80_E2 0.00669438002290
#define GRS80_GM 3.986005e14
#define GRS80_J2 1.08263e-3
#define GRS80_OMEGA 7.292115e-5
#define MGAL_PER_MS2 1.0e5

#define FUNCTIONAL_TAG "V_T_dg"
#define ENGINE_LABEL "CHarm " CHARM_VERSION_TAG " (native, single algorithm)"
#define TIMINGS_HEADER "engine,recurrence,set,nmax,n_points,threads,repeat,wall_s,cpu_s,notes"
#define RESULT_HEADER "id,lat_deg,lon_deg,h_m,V_m2s2,T_m2s2,dg_mgal"

#define CHECK(err)                                   \
    do                                               \
    {                                                \
        if (!charm_err_isempty(err))                 \
            charm_err_handler(err, 1);               \
    } while (0)

/* ------------------------------------------------------------------------------------------------------------ */
/* Timers: wall (monotonic) and process CPU time (all threads)                                                   */
/* ------------------------------------------------------------------------------------------------------------ */
static double wall_now(void)
{
#ifdef _WIN32
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

static double cpu_now(void)
{
#ifdef _WIN32
    FILETIME c, e, k, u;
    ULARGE_INTEGER kk, uu;
    GetProcessTimes(GetCurrentProcess(), &c, &e, &k, &u);
    kk.LowPart = k.dwLowDateTime;
    kk.HighPart = k.dwHighDateTime;
    uu.LowPart = u.dwLowDateTime;
    uu.HighPart = u.dwHighDateTime;
    return 1e-7 * (double)(kk.QuadPart + uu.QuadPart);
#else
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

static void die(const char *msg)
{
    fprintf(stderr, "charm_shs_bench: %s\n", msg);
    exit(1);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Command line                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */
typedef struct
{
    const char *gfc, *points, *out, *timings, *set, *parallel;
    long nmax;
    int threads, repeat;
} Args;

static void usage(const char *why)
{
    fprintf(stderr, "charm_shs_bench: %s\n", why);
    fprintf(stderr, "usage: charm_shs_bench --gfc <file.gfc> --points <csv> --nmax N --out <csv>\n"
                    "                       --timings <csv> --set colorado|latsweep --threads K --repeat R\n"
                    "                       [--parallel points|charm]   (default charm)\n");
    exit(2);
}

static Args parse_args(int argc, char **argv)
{
    Args a;
    memset(&a, 0, sizeof a);
    a.nmax = -1;
    a.threads = 0;
    a.repeat = 1;
    a.parallel = "charm";
    for (int i = 1; i < argc; i++)
    {
        const char *key = argv[i];
        if (i + 1 >= argc)
            usage("missing value after option");
        const char *val = argv[++i];
        if (strcmp(key, "--gfc") == 0)
            a.gfc = val;
        else if (strcmp(key, "--points") == 0)
            a.points = val;
        else if (strcmp(key, "--nmax") == 0)
            a.nmax = strtol(val, NULL, 10);
        else if (strcmp(key, "--out") == 0)
            a.out = val;
        else if (strcmp(key, "--timings") == 0)
            a.timings = val;
        else if (strcmp(key, "--set") == 0)
            a.set = val;
        else if (strcmp(key, "--threads") == 0)
            a.threads = atoi(val);
        else if (strcmp(key, "--repeat") == 0)
            a.repeat = atoi(val);
        else if (strcmp(key, "--recurrence") == 0)
        {
            if (strcmp(val, "native") != 0)
                usage("--recurrence must be native (CHarm has one algorithm)");
        }
        else if (strcmp(key, "--parallel") == 0)
            a.parallel = val;
        else
            usage("unknown option");
    }
    if (!a.gfc || !a.points || !a.out || !a.timings)
        usage("--gfc, --points, --out and --timings are required");
    if (!a.set || (strcmp(a.set, "colorado") != 0 && strcmp(a.set, "latsweep") != 0))
        usage("--set must be colorado or latsweep");
    if (a.nmax < 0)
        usage("--nmax must be a non-negative integer");
    if (a.threads < 1 || a.repeat < 1)
        usage("--threads and --repeat must be >= 1");
    if (strcmp(a.parallel, "points") != 0 && strcmp(a.parallel, "charm") != 0)
        usage("--parallel must be points or charm");
    return a;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Header of the .gfc (same keys and checks as code/python/gfc_reader.py)                                       */
/* ------------------------------------------------------------------------------------------------------------ */
typedef struct
{
    char modelname[256];
    char norm[64];
    char tide_system[64];
    double gm, radius;
    long max_degree;
} GfcHeader;

static char *trim(char *s)
{
    while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n')
        s++;
    char *e = s + strlen(s);
    while (e > s && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\r' || e[-1] == '\n'))
        *--e = '\0';
    return s;
}

static void lower_inplace(char *s)
{
    for (; *s; s++)
        *s = (char)((*s >= 'A' && *s <= 'Z') ? *s + 32 : *s);
}

/* Fortran exponents (D) are accepted, as in gfc_reader.py. */
static double parse_gfc_double(const char *text, const char *path)
{
    char buf[128];
    size_t n = strlen(text);
    if (n >= sizeof buf)
        die("numeric header value too long");
    for (size_t i = 0; i <= n; i++)
        buf[i] = (text[i] == 'D' || text[i] == 'd') ? 'E' : text[i];
    char *end = NULL;
    double v = strtod(buf, &end);
    if (end == buf)
    {
        fprintf(stderr, "charm_shs_bench: %s: bad numeric header value '%s'\n", path, text);
        exit(1);
    }
    return v;
}

static GfcHeader read_gfc_header(const char *path)
{
    FILE *fh = fopen(path, "r");
    if (!fh)
    {
        fprintf(stderr, "charm_shs_bench: cannot open %s\n", path);
        exit(1);
    }
    char gm_text[128] = "", r_text[128] = "", deg_text[128] = "", norm[128] = "", tide[128] = "", model[256] = "";
    int have_gm = 0, have_gm_alt = 0, found_end = 0;
    char line[4096];
    while (fgets(line, sizeof line, fh))
    {
        char *t = trim(line);
        if (*t == '\0')
            continue;
        if (strncmp(t, "end_of_head", 11) == 0)
        {
            found_end = 1;
            break;
        }
        char *sp = t;
        while (*sp && *sp != ' ' && *sp != '\t')
            sp++;
        if (*sp == '\0')
            continue; /* key without value */
        *sp = '\0';
        char *val = trim(sp + 1);
        if (strcmp(t, "earth_gravity_constant") == 0)
        {
            strncpy(gm_text, val, sizeof gm_text - 1);
            have_gm = 1;
        }
        else if (strcmp(t, "gravity_constant") == 0)
        {
            strncpy(gm_text, val, sizeof gm_text - 1);
            have_gm_alt = 1;
        }
        else if (strcmp(t, "radius") == 0)
            strncpy(r_text, val, sizeof r_text - 1);
        else if (strcmp(t, "max_degree") == 0)
            strncpy(deg_text, val, sizeof deg_text - 1);
        else if (strcmp(t, "norm") == 0)
            strncpy(norm, val, sizeof norm - 1);
        else if (strcmp(t, "tide_system") == 0)
            strncpy(tide, val, sizeof tide - 1);
        else if (strcmp(t, "modelname") == 0)
            strncpy(model, val, sizeof model - 1);
    }
    fclose(fh);
    if (!found_end)
    {
        fprintf(stderr, "charm_shs_bench: %s: no end_of_head line found\n", path);
        exit(1);
    }
    if (!have_gm && !have_gm_alt)
        die("header: earth_gravity_constant (or gravity_constant) missing");
    if (!r_text[0] || !deg_text[0] || !norm[0] || !tide[0])
        die("header: radius, max_degree, norm or tide_system missing");

    GfcHeader h;
    memset(&h, 0, sizeof h);
    strncpy(h.modelname, model[0] ? model : "unknown", sizeof h.modelname - 1);
    h.gm = parse_gfc_double(gm_text, path);
    h.radius = parse_gfc_double(r_text, path);
    h.max_degree = strtol(deg_text, NULL, 10);
    strncpy(h.norm, norm, sizeof h.norm - 1);
    strncpy(h.tide_system, tide, sizeof h.tide_system - 1);
    lower_inplace(h.norm);
    lower_inplace(h.tide_system);
    if (strcmp(h.norm, "fully_normalized") != 0)
    {
        fprintf(stderr, "charm_shs_bench: %s: norm is '%s', expected fully_normalized\n", path, h.norm);
        exit(1);
    }
    if (strcmp(h.tide_system, "tide_free") != 0)
    {
        fprintf(stderr, "charm_shs_bench: %s: tide_system is '%s', expected tide_free\n", path, h.tide_system);
        exit(1);
    }
    return h;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Points CSV: header id,lat_deg,lon_deg,h_m                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */
typedef struct
{
    long n;
    char **ids;
    double *lat, *lon, *h;
} Points;

static void read_points(const char *path, Points *p)
{
    FILE *fh = fopen(path, "r");
    if (!fh)
    {
        fprintf(stderr, "charm_shs_bench: cannot open %s\n", path);
        exit(1);
    }
    char line[4096];
    int header_ok = 0;
    long cap = 0;
    p->n = 0;
    p->ids = NULL;
    p->lat = p->lon = p->h = NULL;
    while (fgets(line, sizeof line, fh))
    {
        char *t = trim(line);
        if (*t == '\0')
            continue;
        if (!header_ok)
        {
            if (strcmp(t, "id,lat_deg,lon_deg,h_m") != 0)
            {
                fprintf(stderr, "charm_shs_bench: %s: header must be id,lat_deg,lon_deg,h_m\n", path);
                exit(1);
            }
            header_ok = 1;
            continue;
        }
        /* split the row on commas (no quoting in the benchmark inputs) */
        char *f[4];
        int nf = 0;
        char *s = t;
        f[nf++] = s;
        for (char *c = s; *c; c++)
            if (*c == ',')
            {
                *c = '\0';
                if (nf < 4)
                    f[nf] = c + 1;
                nf++;
            }
        if (nf < 4)
        {
            fprintf(stderr, "charm_shs_bench: %s: row %ld has fewer than 4 fields\n", path, p->n + 2);
            exit(1);
        }
        if (p->n == cap)
        {
            cap = cap ? 2 * cap : 512;
            p->ids = realloc(p->ids, (size_t)cap * sizeof *p->ids);
            p->lat = realloc(p->lat, (size_t)cap * sizeof *p->lat);
            p->lon = realloc(p->lon, (size_t)cap * sizeof *p->lon);
            p->h = realloc(p->h, (size_t)cap * sizeof *p->h);
            if (!p->ids || !p->lat || !p->lon || !p->h)
                die("out of memory reading points");
        }
        p->ids[p->n] = strdup(trim(f[0]));
        p->lat[p->n] = strtod(trim(f[1]), NULL);
        p->lon[p->n] = strtod(trim(f[2]), NULL);
        p->h[p->n] = strtod(trim(f[3]), NULL);
        if (fabs(p->lat[p->n]) > 90.0)
        {
            fprintf(stderr, "charm_shs_bench: %s: latitude outside [-90, 90] at row %ld\n", path, p->n + 2);
            exit(1);
        }
        p->n++;
    }
    fclose(fh);
    if (!header_ok || p->n == 0)
    {
        fprintf(stderr, "charm_shs_bench: %s: no points\n", path);
        exit(1);
    }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Geometry and normal field: ports of shs_common.py                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

/* geodetic_to_geocentric: r [m], phi' [rad] (atan2(Z, p)), lambda [rad] */
static void geodetic_to_geocentric(double lat_deg, double lon_deg, double h_m, double *r, double *phi_c, double *lam)
{
    const double deg2rad = M_PI / 180.0;
    double phi = lat_deg * deg2rad;
    *lam = lon_deg * deg2rad;
    double sin_phi = sin(phi);
    double n = GRS80_A / sqrt(1.0 - GRS80_E2 * sin_phi * sin_phi);
    double p = (n + h_m) * cos(phi);
    double z = (n * (1.0 - GRS80_E2) + h_m) * sin_phi;
    *r = hypot(p, z);
    *phi_c = atan2(z, p);
}

/* q(x) of normal_potential_closed() in shs_common.py */
static double q_func(double x, double e_lin)
{
    return 0.5 * ((1.0 + 3.0 * x * x / (e_lin * e_lin)) * atan(e_lin / x) - 3.0 * x / e_lin);
}

/* normal_potential_closed(): gravitational part of the GRS80 normal potential [m^2/s^2] */
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
    double flattening = 0.5 * GRS80_OMEGA * GRS80_OMEGA * GRS80_A * GRS80_A * (q_func(u, e_lin) / q0)
                        * (sin2_beta - 1.0 / 3.0);
    return GRS80_GM / e_lin * atan(e_lin / u) + flattening;
}

/* ellipsoid_zonal(): fully normalised zonal coefficient C^Ell_{2k,0} */
static double ellipsoid_zonal(int k)
{
    if (k == 0)
        return 1.0;
    double c20 = -GRS80_J2 / sqrt(5.0);
    if (k == 1)
        return c20;
    double e2 = GRS80_E2;
    double sign = (k % 2 == 0) ? 1.0 : -1.0;
    return sign * 3.0 * pow(e2, (double)k) / ((2.0 * k + 1.0) * (2.0 * k + 3.0) * sqrt(4.0 * k + 1.0))
           * (1.0 - k - 5.0 * sqrt(5.0) * k * c20 / e2);
}

/* normal_gravity_gamma(): gamma_SH(r, phi') [m/s^2], zonal series up to degree 20 plus the centrifugal term */
static double normal_gravity_gamma(double r, double phi_c)
{
    double t = sin(phi_c);
    double u = cos(phi_c);
    if (fabs(phi_c) == M_PI / 2.0)
    {
        t = (phi_c > 0.0) ? 1.0 : -1.0; /* np.sign(phi_c) */
        u = 0.0;
    }
    double q = GRS80_A / r;

    double p0[21], p1[21];
    p0[0] = 1.0;
    p0[1] = sqrt(3.0) * t;
    p1[0] = 0.0;
    p1[1] = sqrt(3.0) * u;
    for (int n = 2; n <= 20; n++)
    {
        double nn = (double)n;
        p0[n] = sqrt(4.0 * nn * nn - 1.0) / nn * t * p0[n - 1]
                - (nn - 1.0) / nn * sqrt((2.0 * nn + 1.0) / (2.0 * nn - 3.0)) * p0[n - 2];
        double a1 = sqrt((2.0 * nn + 1.0) * (2.0 * nn - 1.0) / ((nn - 1.0) * (nn + 1.0)));
        double b1 = (n >= 3) ? sqrt((2.0 * nn + 1.0) * nn * (nn - 2.0) / ((nn - 1.0) * (nn + 1.0) * (2.0 * nn - 3.0)))
                             : 0.0;
        p1[n] = a1 * t * p1[n - 1] - b1 * p1[n - 2];
    }

    double sr = 0.0, sp = 0.0, qn = 1.0;
    for (int k = 0; k <= 10; k++)
    {
        int n = 2 * k;
        double c = ellipsoid_zonal(k);
        sr += qn * (n + 1) * c * p0[n];
        sp += qn * c * sqrt(n * (n + 1) / 2.0) * p1[n];
        qn = qn * q * q;
    }
    double dvr = -GRS80_GM / (r * r) * sr;
    double dvp = GRS80_GM / r * sp;
    double w2 = GRS80_OMEGA * GRS80_OMEGA;
    double dvcr = w2 * r * u * u;
    double dvcp = -w2 * r * r * u * t;
    double g_x = (dvp + dvcp) / r;
    double g_z = dvr + dvcr;
    return hypot(g_x, g_z);
}

/* gravity_disturbance_mgal(): dg = (|g| - gamma_SH) * 1e5, with the centrifugal gradient added here */
static double gravity_disturbance_mgal(double v_r, double v_lat, double v_lon, double r, double phi_c)
{
    double cos_p = cos(phi_c);
    double sin_p = sin(phi_c);
    double w2 = GRS80_OMEGA * GRS80_OMEGA;
    double vc_r = w2 * r * cos_p * cos_p;
    double vc_lat = -w2 * r * r * cos_p * sin_p;
    double g_x = (v_lat + vc_lat) / r;
    double g_y = -v_lon / (r * cos_p);
    double g_z = v_r + vc_r;
    double g = sqrt(g_x * g_x + g_y * g_y + g_z * g_z);
    return (g - normal_gravity_gamma(r, phi_c)) * MGAL_PER_MS2;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Synthesis                                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/* One point per iteration: the point structure borrows the driver's arrays (owner = 0), so each thread owns its
 * own charm_point and charm_err. The shared charm_shc is read-only during the synthesis. */
static void synthesise_points(const charm_shc *shcs, unsigned long nmax, long n, double *lat, double *lon,
                              double *rr, double *v, double *gp, double *gl, double *gr, int threads, int batch)
{
    if (batch)
    {
        /* --parallel charm: one CHarm call over all points (the way pyharm calls CHarm). CHarm's own OpenMP
         * region, over the orders m, runs with 'threads' threads; the driver adds no parallelism here. */
        omp_set_num_threads(threads);
        charm_err *err = charm_err_init();
        if (!err)
        {
            fprintf(stderr, "charm_shs_bench: charm_err_init failed\n");
            exit(1);
        }
        charm_point *pnt = charm_crd_point_init(CHARM_CRD_POINT_SCATTERED, (size_t)n, (size_t)n, lat, lon, rr);
        if (!pnt)
        {
            fprintf(stderr, "charm_shs_bench: charm_crd_point_init failed\n");
            exit(1);
        }
        charm_shs_point(pnt, shcs, nmax, v, err);
        CHECK(err);
        /* the gradient buffers: charm_shs_point_grad1 swaps the pointers, so keep the originals for freeing */
        double *b0 = malloc((size_t)n * sizeof *b0), *b1 = malloc((size_t)n * sizeof *b1),
               *b2 = malloc((size_t)n * sizeof *b2);
        if (!b0 || !b1 || !b2)
        {
            fprintf(stderr, "charm_shs_bench: out of memory\n");
            exit(1);
        }
        double *fg[3] = { b0, b1, b2 };
        charm_shs_point_grad1(pnt, shcs, nmax, fg, err);
        CHECK(err);
        for (long i = 0; i < n; i++)
        {
            gp[i] = fg[0][i];
            gl[i] = fg[1][i];
            gr[i] = fg[2][i];
        }
        free(b0);
        free(b1);
        free(b2);
        charm_crd_point_free(pnt);
        charm_err_free(err);
        return;
    }

#pragma omp parallel num_threads(threads) default(shared)
    {
        /* CHarm's internal region (over orders m) must run with one thread inside each driver thread. With
         * num_threads(1) the outer region is not active, so max-active-levels alone does not stop the inner
         * region from taking the default team size; setting the nested thread count per thread does. */
        omp_set_num_threads(1);

#pragma omp for schedule(dynamic, 4)
        for (long i = 0; i < n; i++)
        {
            charm_err *err = charm_err_init();
            if (!err)
            {
                fprintf(stderr, "charm_shs_bench: charm_err_init failed\n");
                exit(1);
            }
            charm_point *pnt = charm_crd_point_init(CHARM_CRD_POINT_SCATTERED, 1, 1, &lat[i], &lon[i], &rr[i]);
            if (!pnt)
            {
                fprintf(stderr, "charm_shs_bench: charm_crd_point_init failed at point %ld\n", i);
                exit(1);
            }
            charm_shs_point(pnt, shcs, nmax, &v[i], err);
            CHECK(err);
            /* charm_shs_point_grad1 swaps the pointers of fg on return: after the call fg[0], fg[1], fg[2] hold
             * the LNOF x (latitude, 1/r dV/dphi), y (longitude, -1/(r cos phi) dV/dlambda) and z (radial, dV/dr)
             * elements, so the values are read through fg, not through the buffers passed in (verified against
             * finite differences of V and against pyharm, see README.md). */
            double gx = 0.0, gy = 0.0, gz = 0.0;
            double *fg[3] = { &gx, &gy, &gz };
            charm_shs_point_grad1(pnt, shcs, nmax, fg, err);
            CHECK(err);
            gp[i] = *fg[0];
            gl[i] = *fg[1];
            gr[i] = *fg[2];
            charm_crd_point_free(pnt);
            charm_err_free(err);
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Output                                                                                                        */
/* ------------------------------------------------------------------------------------------------------------ */

/* Commas would break the CSV timings row: replace them by semicolons in free text. */
static void sanitize(char *s)
{
    for (; *s; s++)
        if (*s == ',' || *s == '"' || *s == '\n' || *s == '\r')
            *s = ';';
}

static void write_result_csv(const char *path, const Points *pts, const double *v, const double *t,
                             const double *dg)
{
    FILE *fh = fopen(path, "w");
    if (!fh)
    {
        fprintf(stderr, "charm_shs_bench: cannot write %s\n", path);
        exit(1);
    }
    fprintf(fh, RESULT_HEADER "\n");
    for (long i = 0; i < pts->n; i++)
        fprintf(fh, "%s,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n", pts->ids[i], pts->lat[i], pts->lon[i], pts->h[i],
                v[i], t[i], dg[i]);
    fclose(fh);
}

static void make_md_path(const char *out, char *md, size_t cap)
{
    size_t n = strlen(out);
    if (n > 4 && strcmp(out + n - 4, ".csv") == 0)
        n -= 4;
    if (n + 4 >= cap)
        die("output path too long");
    memcpy(md, out, n);
    strcpy(md + n, ".md");
}

static const char *basename_of(const char *path)
{
    const char *b = path, *p;
    for (p = path; *p; p++)
        if (*p == '/' || *p == '\\')
            b = p + 1;
    return b;
}

static void write_result_md(const char *md_path, const char *out, const char *command, const Args *a,
                            const GfcHeader *hdr, double gm_model, double r_model, double load_s, long n_points,
                            const char *compiler, const char *threads_note)
{
    FILE *fh = fopen(md_path, "w");
    if (!fh)
    {
        fprintf(stderr, "charm_shs_bench: cannot write %s\n", md_path);
        exit(1);
    }
    time_t now = time(NULL);
    char stamp[64];
    strftime(stamp, sizeof stamp, "%Y-%m-%dT%H:%M:%S", localtime(&now));
    fprintf(fh, "# %s\n\n", basename_of(out));
    fprintf(fh, "Engine: **%s**. Generated by the runner script on %s. V, T and dg are computed in the same run.\n\n",
            ENGINE_LABEL, stamp);
    fprintf(fh, "## Columns\n\n| column | unit | meaning |\n|---|---|---|\n");
    fprintf(fh, "| id | - | point identifier from the input CSV |\n");
    fprintf(fh, "| lat_deg | deg | geodetic latitude on GRS80 (input, echoed) |\n");
    fprintf(fh, "| lon_deg | deg | longitude (input, echoed) |\n");
    fprintf(fh, "| h_m | m | ellipsoidal height on GRS80 (input, echoed) |\n");
    fprintf(fh, "| V_m2s2 | m^2/s^2 | gravitational potential V of the model at the point (no normal potential "
                "subtracted; degree 0 included, C00 as in the file) |\n");
    fprintf(fh, "| T_m2s2 | m^2/s^2 | disturbing potential T = V - U at the point (V: model, U: GRS80 normal "
                "potential, closed form) |\n");
    fprintf(fh, "| dg_mgal | mGal | gravity disturbance dg = \\|g\\| - gamma_SH at the point (g includes the "
                "centrifugal term; gamma_SH: project series) |\n\n");
    fprintf(fh, "## Command and parameters\n\n```\n%s\n```\n\n", command);
    fprintf(fh, "- nmax: `%ld`\n", a->nmax);
    fprintf(fh, "- set: `%s`\n", a->set);
    fprintf(fh, "- repeat: `%d`\n", a->repeat);
    fprintf(fh, "- threads: `%d`\n", a->threads);
    fprintf(fh, "- normalization: `4pi (fully normalised)`\n");
    fprintf(fh, "- csphase: `none (GFC convention)`\n");
    fprintf(fh, "- functionals: `V, T, dg from the same synthesis pass`\n");
    fprintf(fh, "- gravity disturbance: `|g| - gamma_SH with centrifugal term, omega=7.292115e-5`\n");
    fprintf(fh, "- parallelisation: `%s`\n\n", threads_note);
    fprintf(fh, "## Model\n\n");
    fprintf(fh, "- file: `%s`\n", basename_of(a->gfc));
    fprintf(fh, "- modelname: `%s`\n", hdr->modelname);
    fprintf(fh, "- GM [m^3/s^2]: `%.10e`\n", gm_model);
    fprintf(fh, "- R [m]: `%.6f`\n", r_model);
    fprintf(fh, "- tide_system: `%s`\n", hdr->tide_system);
    fprintf(fh, "- load_s: `%.3f`\n\n", load_s);
    fprintf(fh, "Points: %ld.\n\n", n_points);
    fprintf(fh, "## Notes\n\n");
    fprintf(fh, "- V = charm_shs_point: the model potential itself, no normal potential subtracted, degree 0 "
                "included (C00 as in the file).\n");
    fprintf(fh, "- T = V - U, U = GRS80 level-ellipsoid gravitational potential (closed form), with the V of this "
                "run.\n");
    fprintf(fh, "- Functionals evaluated at the actual point (geodetic h), not on the ellipsoid.\n");
    fprintf(fh, "- Gradient from charm_shs_point_grad1 (LNOF, latitude component divided by r); converted back to "
                "partial derivatives of V as in run_pyharm.py.\n");
    fprintf(fh, "- CHarm %s; compiler %s; flags %s.\n", CHARM_VERSION_TAG, compiler, CHARM_DRIVER_FLAGS);
    fclose(fh);
}

static void append_timing(const char *path, const Args *a, long n_points, int repeat, double wall, double cpu,
                          const char *notes)
{
    FILE *probe = fopen(path, "r");
    int is_new = 1;
    if (probe)
    {
        fseek(probe, 0, SEEK_END);
        is_new = (ftell(probe) == 0);
        fclose(probe);
    }
    FILE *fh = fopen(path, "a");
    if (!fh)
    {
        fprintf(stderr, "charm_shs_bench: cannot append to %s\n", path);
        exit(1);
    }
    if (is_new)
        fprintf(fh, TIMINGS_HEADER "\n");
    fprintf(fh, "charm,native,%s,%ld,%ld,%d,%d,%.6f,%.6f,%s\n", a->set, a->nmax, n_points, a->threads, repeat, wall,
            cpu, notes);
    fclose(fh);
}

static void compiler_string(char *buf, size_t cap)
{
#if defined(__VERSION__)
    snprintf(buf, cap, "%s", __VERSION__);
#else
    snprintf(buf, cap, "unknown");
#endif
    sanitize(buf);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* main                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    Args a = parse_args(argc, argv);
    GfcHeader hdr = read_gfc_header(a.gfc);
    if (a.nmax > hdr.max_degree)
    {
        fprintf(stderr, "charm_shs_bench: nmax=%ld exceeds max_degree=%ld of the model\n", a.nmax, hdr.max_degree);
        return 1;
    }

    /* Nested OpenMP off: CHarm's internal region (over orders m) runs with one thread inside the point loop. */
    omp_set_dynamic(0);
    omp_set_max_active_levels(1);

    /* ---- model load (timed separately, reported in notes) ---- */
    double t0 = wall_now();
    charm_err *err = charm_err_init();
    if (!err)
        die("charm_err_init failed");
    charm_shc *shcs = charm_shc_calloc((unsigned long)a.nmax, 1.0, 1.0);
    if (!shcs)
        die("charm_shc_calloc failed");
    charm_shc_read_gfc(a.gfc, (unsigned long)a.nmax, NULL, shcs, err);
    CHECK(err);
    double load_s = wall_now() - t0;
    charm_err_free(err);

    double gm_model = shcs->mu, r_model = shcs->r;
    if (fabs(gm_model - hdr.gm) > 1e-9 * hdr.gm || fabs(r_model - hdr.radius) > 1e-6)
    {
        fprintf(stderr, "charm_shs_bench: GM/R mismatch: CHarm (%.10e, %.6f) vs header (%.10e, %.6f)\n", gm_model,
                r_model, hdr.gm, hdr.radius);
        return 1;
    }

    /* ---- points and geometry (identical to shs_common.geodetic_to_geocentric) ---- */
    Points pts;
    read_points(a.points, &pts);
    long n = pts.n;
    double *rr = malloc((size_t)n * sizeof *rr), *phic = malloc((size_t)n * sizeof *phic),
           *lam = malloc((size_t)n * sizeof *lam);
    double *v = malloc((size_t)n * sizeof *v), *gp = malloc((size_t)n * sizeof *gp),
           *gl = malloc((size_t)n * sizeof *gl), *gr = malloc((size_t)n * sizeof *gr);
    double *tv = malloc((size_t)n * sizeof *tv), *dg = malloc((size_t)n * sizeof *dg);
    if (!rr || !phic || !lam || !v || !gp || !gl || !gr || !tv || !dg)
        die("out of memory");
    for (long i = 0; i < n; i++)
        geodetic_to_geocentric(pts.lat[i], pts.lon[i], pts.h[i], &rr[i], &phic[i], &lam[i]);

    /* ---- timed repetitions: synthesis of V and the LNOF gradient at every point ---- */
    char compiler[256];
    compiler_string(compiler, sizeof compiler);
    const int batch = strcmp(a.parallel, "charm") == 0;
    const char *threads_note =
        batch ? "parallel=charm: one CHarm call over all points; CHarm internal OpenMP over orders m with K threads; "
                "the driver adds no parallelism"
              : "parallel=points: OpenMP over points in the driver (K threads, one CHarm call per point); CHarm "
                "internal OpenMP over orders m kept at one thread (nested level 1)";
    const char *call_note = batch ? "charm_shs_point+charm_shs_point_grad1 one call for all points"
                                  : "charm_shs_point+charm_shs_point_grad1 one point per call";
    char notes[1024];
    for (int rep = 1; rep <= a.repeat; rep++)
    {
        double w0 = wall_now();
        double c0 = cpu_now();
        synthesise_points(shcs, (unsigned long)a.nmax, n, phic, lam, rr, v, gp, gl, gr, a.threads, batch);
        double cpu = cpu_now() - c0;
        double wall = wall_now() - w0;
        snprintf(notes, sizeof notes,
                 "load_s=%.3f; %s; %s; functional=%s; charm=%s; compiler=%s; flags=%s; platform=%s",
                 load_s, call_note, threads_note, FUNCTIONAL_TAG, CHARM_VERSION_TAG, compiler, CHARM_DRIVER_FLAGS,
#ifdef _WIN32
                 "Windows"
#else
                 "POSIX"
#endif
        );
        sanitize(notes);
        append_timing(a.timings, &a, n, rep, wall, cpu, notes);
        printf("repeat %d/%d: wall=%.3fs cpu=%.3fs\n", rep, a.repeat, wall, cpu);
        fflush(stdout);
    }

    /* ---- functionals from the last repetition (same formulas as run_pyharm.py) ---- */
    for (long i = 0; i < n; i++)
    {
        double v_r = gr[i];
        double v_lat = gp[i] * rr[i];                           /* dV/dphi' = r f[P] */
        double v_lon = -gl[i] * rr[i] * cos(phic[i]);           /* dV/dlambda = -r cos(phi') f[L] */
        tv[i] = v[i] - normal_potential_closed(rr[i], phic[i]); /* T = V - U */
        dg[i] = gravity_disturbance_mgal(v_r, v_lat, v_lon, rr[i], phic[i]);
    }

    write_result_csv(a.out, &pts, v, tv, dg);
    char md_path[4096];
    make_md_path(a.out, md_path, sizeof md_path);
    char command[8192] = "";
    for (int i = 0; i < argc; i++)
    {
        strncat(command, i ? " " : "", sizeof command - strlen(command) - 1);
        strncat(command, argv[i], sizeof command - strlen(command) - 1);
    }
    write_result_md(md_path, a.out, command, &a, &hdr, gm_model, r_model, load_s, n, compiler, threads_note);

    printf("first rows:\n");
    for (long i = 0; i < n && i < 3; i++)
        printf("%s,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n", pts.ids[i], pts.lat[i], pts.lon[i], pts.h[i], v[i], tv[i],
               dg[i]);
    printf("wrote %s and %s\n", a.out, md_path);

    charm_shc_free(shcs);
    return 0;
}
