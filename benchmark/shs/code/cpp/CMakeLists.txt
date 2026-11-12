# Driver of the gsh engine for the SHS benchmark (definitions: benchmark/shs/results/csv.md).
add_executable(gsh_shs_bench gsh_shs_bench.cpp)
target_link_libraries(gsh_shs_bench PRIVATE gsh)
