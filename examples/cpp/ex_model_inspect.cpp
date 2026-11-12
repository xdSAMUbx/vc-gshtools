// ex_model_inspect.cpp: reads an ICGEM .gfc model (model.hpp) and prints its header, degree variances and low-degree terms.
//
// Shows: load_gfc() and the fields of gsh::Model (GM, R, degrees, norm, tide system), the packed coefficient accessors
// Model::c() and Model::s(), and degree_variances() (errors.hpp).
//
// Usage: ex_model_inspect <model.gfc> [nmax]
//   nmax: optional maximum degree to read (>= 2). Without it the whole file is read (EGM2008: degree 2190).
//
// Output: header fields; the degree variances sigma_n^2 = sum_m (C_nm^2 + S_nm^2) for a few degrees, and the
// degree amplitude sqrt(sigma_n^2), which is the typical size of the coefficients of degree n (dimensionless,
// fully normalised); and C20, C21 and S21 in the same normalisation.

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include "errors.hpp"
#include "model.hpp"

namespace
{
  [[noreturn]] void usage(const char *why)
  {
    std::fprintf(stderr, "ex_model_inspect: %s\nusage: ex_model_inspect <model.gfc> [nmax]\n", why);
    std::exit(2);
  }

  int parse_int(const char *s)
  {
    char *end = nullptr;
    const long v = std::strtol(s, &end, 10);
    if (end == s || *end != '\0')
      usage("nmax must be an integer");
    return static_cast<int>(v);
  }
} // namespace

int main(int argc, char **argv)
{
  try
  {
    if (argc < 2 || argc > 3)
      usage("expected 1 or 2 arguments");
    const int nmax = argc > 2 ? parse_int(argv[2]) : -1;
    if (argc > 2 && nmax < 2)
      usage("nmax must be at least 2");

    const gsh::Model model = gsh::load_gfc(argv[1], nmax);
    if (model.nmax < 2)
      throw std::runtime_error("the model has no coefficients of degree 2 or higher");

    std::printf("Model: %s\n", argv[1]);
    std::printf("  name:        %s\n", model.name.c_str());
    std::printf("  GM:          %.10e m^3/s^2\n", model.GM);
    std::printf("  R:           %.8f m\n", model.R);
    std::printf("  file degree: %d, degrees read: %d\n", model.file_nmax, model.nmax);
    std::printf("  norm:        %s\n", model.norm.c_str());
    std::printf("  tide system: %s\n", model.tide_system.c_str());

    const std::vector<double> var = gsh::degree_variances(model);
    std::printf("Degree variances sigma_n^2 and amplitudes sqrt(sigma_n^2):\n");
    for (int n : {2, 10, 100, 1000, 2190})
    {
      if (n > model.nmax)
        continue;
      std::printf("  n = %4d: sigma^2 = %.6e, sqrt = %.6e\n", n, var[static_cast<std::size_t>(n)],
                  std::sqrt(var[static_cast<std::size_t>(n)]));
    }

    std::printf("Low-degree coefficients (fully normalised):\n");
    std::printf("  C20 = %.12e\n", model.c(2, 0));
    std::printf("  C21 = %.12e\n", model.c(2, 1));
    std::printf("  S21 = %.12e\n", model.s(2, 1));
    return 0;
  }
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "ex_model_inspect: %s\n", e.what());
    return 1;
  }
}
