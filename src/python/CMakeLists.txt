# src/python/CMakeLists.txt: pybind11 extension gshtools._core, linked to the C++ library target gsh.
# Enabled from the root with -DGSH_BUILD_PYTHON=ON. Release is always optimised (root guard for /O2 applies here).

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core
  src/module.cpp
  src/common.hpp
  src/model_points.cpp
  src/geo.cpp
  src/legendre.cpp
  src/functionals.cpp
  src/synthesis.cpp
  src/errors.cpp
)
target_link_libraries(_core PRIVATE gsh)
set_target_properties(_core PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gshtools"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gshtools")

# The extension needs no third-party DLL (the FFT is the vendored, header-only pocketfft). The bundler still takes a
# list of extra search directories, which is empty.
set(_gsh_extra_dirs "")

# In-tree layout: build/<dir>/python/gshtools/ holds the extension, __init__.py and the DLLs it needs, so that
# PYTHONPATH=<build dir>/python makes `import gshtools` work without installing.
add_custom_command(TARGET _core POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/gshtools/__init__.py" "${CMAKE_CURRENT_BINARY_DIR}/gshtools/__init__.py"
  COMMAND ${CMAKE_COMMAND}
          "-DGSH_PYD=$<TARGET_FILE:_core>" "-DGSH_DEST=$<TARGET_FILE_DIR:_core>"
          "-DGSH_EXTRA_DIRS=${_gsh_extra_dirs}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bundle_runtime.cmake"
  VERBATIM)

# Installed layout (pip install . or cmake --install): gshtools/ with __init__.py, the extension and its runtime DLLs.
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/gshtools/__init__.py" DESTINATION gshtools)
install(TARGETS _core LIBRARY DESTINATION gshtools RUNTIME DESTINATION gshtools)
install(CODE "set(GSH_EXTRA_DIRS [[${_gsh_extra_dirs}]])\ninclude([[${CMAKE_CURRENT_SOURCE_DIR}/cmake/bundle_runtime.cmake]])")
