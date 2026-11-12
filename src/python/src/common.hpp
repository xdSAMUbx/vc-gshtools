// common.hpp: conversion helpers shared by the binding sources. The wrapper only converts types and calls the
// C++ library; no numerical logic lives here.
//
// Memory policy:
//   - Inputs: the C++ API takes owning containers (Points, std::vector), so every array input is copied into one.
//     Arrays that are already float64 and C-contiguous are not converted before that copy.
//   - Outputs: zero-copy. The C++ container is moved to the heap and NumPy takes over its storage through a
//     capsule that frees it when the array is collected.
//   - Read-only views (for example Points.lat) reference the owning Python object and copy nothing.
#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include "core.hpp"

namespace gshpy
{
  namespace py = pybind11;
  using ssize = py::ssize_t;

  // Float64 input: C-contiguous; other dtypes and layouts are converted (a copy).
  using DArr = py::array_t<double, py::array::c_style | py::array::forcecast>;

  // Copies a one-dimensional input into an owning container V (AlignedVec<double> or std::vector<double>).
  template <class V>
  V copy_in(const DArr &a, const char *name)
  {
    if (a.ndim() != 1)
      throw py::value_error(std::string(name) + " must be one-dimensional (got ndim=" + std::to_string(a.ndim()) + ")");
    const double *p = a.data();
    return V(p, p + a.size());
  }

  // Throws ValueError when two inputs that must have the same length do not.
  inline void require_same(ssize a, ssize b, const char *what)
  {
    if (a != b)
      throw py::value_error(std::string(what) + ": the arrays must have the same length");
  }

  // Moves an owning container to NumPy without copying its elements.
  template <class V>
  py::array to_numpy(V v)
  {
    using T = typename V::value_type;
    V *holder = new V(std::move(v));
    py::capsule owner(holder, [](void *p) { delete static_cast<V *>(p); });
    return py::array_t<T>(std::vector<ssize>{static_cast<ssize>(holder->size())},
                          std::vector<ssize>{static_cast<ssize>(sizeof(T))}, holder->data(), owner);
  }

  // Read-only NumPy view of a member container. `owner` (the Python object that holds the container) is kept alive
  // by the view.
  template <class V>
  py::array view_of(const V &v, py::handle owner)
  {
    using T = typename V::value_type;
    py::array_t<T> a(std::vector<ssize>{static_cast<ssize>(v.size())}, std::vector<ssize>{static_cast<ssize>(sizeof(T))},
                     v.data(), owner);
    a.attr("setflags")(py::arg("write") = false);
    return a;
  }

  // Bind functions, one per source file. Each takes the module (or submodule) that owns the names it registers.
  void bind_model_points(py::module_ &m);  // Model, load_gfc, Points, SynthOptions, Isa, read_table
  void bind_ellipsoid(py::module_ &m);     // Ellipsoid, disturbing_model
  void bind_coords(py::module_ &m);        // submodule coords
  void bind_tides(py::module_ &m);         // submodule tides (TideSystem included)
  void bind_legendre(py::module_ &m);      // submodule legendre
  void bind_functionals(py::module_ &m);   // Functional, FunctionalRequest, compute_functionals, height sources
  void bind_synthesis(py::module_ &m);     // Prepared, RawFields, synth_points, Grid, synth_grid
  void bind_errors(py::module_ &m);        // errors and error_report
}
