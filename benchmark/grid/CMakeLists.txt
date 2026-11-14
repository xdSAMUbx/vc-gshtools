# Drivers of the gsh side of the grid benchmark (definitions: benchmark/grid/README.md).
add_executable(gsh_grid_bench gsh_grid_bench.cpp)
target_link_libraries(gsh_grid_bench PRIVATE gsh)
if (WIN32)
  # GetProcessMemoryInfo (peak working set) lives in psapi.
  target_link_libraries(gsh_grid_bench PRIVATE psapi)
endif()

add_executable(gsh_grid_check gsh_grid_check.cpp)
target_link_libraries(gsh_grid_check PRIVATE gsh)
