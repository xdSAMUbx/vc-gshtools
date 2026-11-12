# Library target "gsh": utilities, Legendre recurrences with SIMD kernels, and spherical-harmonic synthesis.
add_library(gsh STATIC
  # gsh_utils: constants, model reader, ellipsoid, height, tables, tide systems, synthesis types, coordinates
  gsh_utils/core.cpp
  gsh_utils/model.cpp
  gsh_utils/ellipsoid.cpp
  gsh_utils/height.cpp
  gsh_utils/table_io.cpp
  gsh_utils/tides.cpp
  gsh_utils/tide_system.cpp
  gsh_utils/synth_types.cpp
  gsh_utils/coordinates.cpp
  # gsh_legendre: Legendre functions and internal SIMD kernels
  gsh_legendre/legendre.cpp
  gsh_legendre/kernels/kernels_scalar.cpp
  gsh_legendre/kernels/kernels_avx2.cpp
  gsh_legendre/kernels/kernels_avx512.cpp
  # gsh_expansion: synthesis at points and grids, functionals, height source from a model, reference oracle
  gsh_expansion/shs.cpp
  gsh_expansion/functionals.cpp
  gsh_expansion/grid.cpp
  gsh_expansion/height_model.cpp
  gsh_expansion/reference.cpp
  gsh_expansion/errors.cpp
  gsh_expansion/error_report.cpp
)
target_include_directories(gsh
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gsh_utils
         ${CMAKE_CURRENT_SOURCE_DIR}/gsh_legendre
         ${CMAKE_CURRENT_SOURCE_DIR}/gsh_expansion
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/gsh_legendre/kernels)
target_link_libraries(gsh PUBLIC OpenMP::OpenMP_CXX Eigen3::Eigen PRIVATE FFTW3::fftw3)
target_compile_definitions(gsh PUBLIC _CRT_SECURE_NO_WARNINGS NOMINMAX)

# One translation unit per ISA, so that runtime dispatch is safe.
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(gsh PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/utf-8 /W3 /EHsc>)
  target_compile_options(gsh PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:/O2 /Oi /Ob3 /fp:precise>)
  set_source_files_properties(gsh_legendre/kernels/kernels_avx2.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
  set_source_files_properties(gsh_legendre/kernels/kernels_avx512.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX512")
else() # clang-cl
  target_compile_options(gsh PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/utf-8 /EHsc>)
  target_compile_options(gsh PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:/O2 /Oi>)
  set_source_files_properties(gsh_legendre/kernels/kernels_avx2.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
  set_source_files_properties(gsh_legendre/kernels/kernels_avx512.cpp PROPERTIES COMPILE_OPTIONS "-mavx512f;-mavx512dq;-mavx512vl;-mavx512bw;-mfma;-mprefer-vector-width=512")
endif()
