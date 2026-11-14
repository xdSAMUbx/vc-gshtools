# GRAV-D flight at degree 5540, all 45 functionals (definitions: benchmark/aero5540/README.md).
add_executable(gsh_aero5540 gsh_aero5540.cpp)
target_link_libraries(gsh_aero5540 PRIVATE gsh)
if (WIN32)
  target_link_libraries(gsh_aero5540 PRIVATE psapi)
endif()
