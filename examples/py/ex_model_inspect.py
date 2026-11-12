"""ex_model_inspect.py: reads an ICGEM .gfc model and prints its header, degree variances and low-degree terms (gshtools).

Python equivalent of examples/cpp/ex_model_inspect.cpp. Same arguments, same printed quantities.

Shows: load_gfc() and the fields of Model (GM, R, degrees, norm, tide system), the accessors Model.c() and Model.s(),
and degree_variances().

Usage: python examples/py/ex_model_inspect.py <model.gfc> [nmax]
  nmax: optional maximum degree to read (>= 2). Without it the whole file is read (EGM2008: degree 2190).

Output: header fields; the degree variances sigma_n^2 = sum_m (C_nm^2 + S_nm^2) for a few degrees, and the degree
amplitude sqrt(sigma_n^2), which is the typical size of the coefficients of degree n (dimensionless, fully normalised);
and C20, C21 and S21 in the same normalisation.
"""

import argparse
import math
import sys

import gshtools as gsh


def main():
    parser = argparse.ArgumentParser(description="Header, degree variances and low-degree terms of a .gfc model.")
    parser.add_argument("model", help="ICGEM .gfc model")
    parser.add_argument("nmax", nargs="?", type=int, default=None, help="optional maximum degree to read, >= 2")
    args = parser.parse_args()
    if args.nmax is not None and args.nmax < 2:
        parser.error("nmax must be at least 2")
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    model = gsh.load_gfc(args.model, -1 if args.nmax is None else args.nmax)
    if model.nmax < 2:
        raise RuntimeError("the model has no coefficients of degree 2 or higher")

    print(f"Model: {args.model}")
    print(f"  name:        {model.name}")
    print(f"  GM:          {model.GM:.10e} m^3/s^2")
    print(f"  R:           {model.R:.8f} m")
    print(f"  file degree: {model.file_nmax}, degrees read: {model.nmax}")
    print(f"  norm:        {model.norm}")
    print(f"  tide system: {model.tide_system}")

    var = gsh.degree_variances(model)
    print("Degree variances sigma_n^2 and amplitudes sqrt(sigma_n^2):")
    for n in (2, 10, 100, 1000, 2190):
        if n > model.nmax:
            continue
        print(f"  n = {n:4d}: sigma^2 = {var[n]:.6e}, sqrt = {math.sqrt(var[n]):.6e}")

    print("Low-degree coefficients (fully normalised):")
    print(f"  C20 = {model.c(2, 0):.12e}")
    print(f"  C21 = {model.c(2, 1):.12e}")
    print(f"  S21 = {model.s(2, 1):.12e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
