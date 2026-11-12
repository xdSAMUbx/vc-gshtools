# Driver of the gsh side of the meridional benchmark against GrafLab (definitions: benchmark/meridian/README.md).
add_executable(gsh_meridian gsh_meridian.cpp)
target_link_libraries(gsh_meridian PRIVATE gsh)
