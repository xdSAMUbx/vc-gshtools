# C++ reference of the Python-overhead benchmark (gsh_overhead_ref.cpp, definitions in README.md).
# Added to the build by the root CMakeLists.txt, inside the GSH_BUILD_TOOLS block, with:
#   add_subdirectory(benchmark/python_overhead)
add_executable(gsh_overhead_ref gsh_overhead_ref.cpp)
target_link_libraries(gsh_overhead_ref PRIVATE gsh)
