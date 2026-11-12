// Legendre MSFC en un punto y tiempo: ex_puntual [lat_grados] [nmax]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "legendre.hpp"

int main(int argc, char **argv)
{
  const double lat = argc > 1 ? std::atof(argv[1]) : 10.0;
  const int nmax = argc > 2 ? std::atoi(argv[2]) : 2699;
  const auto coefs = gsh::legendre::make_coefs(nmax);
  std::vector<double> P(gsh::tri_size(nmax));
  const double phi = gsh::deg2rad(lat);
  const auto t0 = std::chrono::steady_clock::now();
  gsh::legendre::compute(std::sin(phi), std::cos(phi), coefs, P);
  const auto t1 = std::chrono::steady_clock::now();
  std::printf("nmax=%d lat=%.3f: %.1f ms, P(nmax,0)=%.15e\n", nmax, lat,
              std::chrono::duration<double, std::milli>(t1 - t0).count(),
              P[gsh::idx_tri(nmax, 0, nmax)]);
}
