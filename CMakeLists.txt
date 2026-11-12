cmake_minimum_required(VERSION 3.24)
project(vc-gshtools VERSION 2.0.0 LANGUAGES CXX
        DESCRIPTION "Vibe-Coded Geodetic Spherical Harmonic Tools: spherical-harmonic synthesis of ultra-high-degree gravity-field models (C++20, SIMD, OpenMP)")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

option(BUILD_TESTING "Build the unit tests (GoogleTest)" OFF)
option(GSH_BUILD_TOOLS "Build the point example and the SHS benchmark driver" ON)

find_package(OpenMP REQUIRED)
find_package(FFTW3 CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# Release must always be optimised, whatever the cache holds. A stale or incomplete CMake cache can leave
# CMAKE_CXX_FLAGS_RELEASE empty, and a clang-cl build without /O2 runs about 40 times slower.
if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
  add_compile_options($<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:/O2>)
endif()

add_subdirectory(src/cpp)

if (GSH_BUILD_TOOLS)
  add_subdirectory(examples)
  add_subdirectory(benchmark/shs/code/cpp)
  add_subdirectory(benchmark/meridian)
endif()

if (BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(benchmark/alf)
endif()

message(STATUS "vc-gshtools: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} | Tests=${BUILD_TESTING} | Tools=${GSH_BUILD_TOOLS}")
