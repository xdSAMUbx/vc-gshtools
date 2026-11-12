# Identity (sum-rule) validation of the Legendre functions and derivatives along the meridian.
# Run: ctest -R AlfIdentities -V   (or build/<preset>/benchmark/alf/gsh_alf_identities)
find_package(GTest CONFIG REQUIRED)
add_executable(gsh_alf_identities test_alf_identities.cpp)
target_link_libraries(gsh_alf_identities PRIVATE gsh GTest::gtest GTest::gtest_main)
target_compile_definitions(gsh_alf_identities PRIVATE GSH_ALF_OUT="${CMAKE_CURRENT_SOURCE_DIR}/results")
include(GoogleTest)
gtest_discover_tests(gsh_alf_identities DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 1800)

# Differential test: production SIMD kernel against the sequential oracle, V along the meridian (EGM2008, 2190).
add_executable(gsh_kernel_vs_oracle test_kernel_vs_oracle.cpp)
target_link_libraries(gsh_kernel_vs_oracle PRIVATE gsh GTest::gtest GTest::gtest_main)
target_compile_definitions(gsh_kernel_vs_oracle PRIVATE
  GSH_ALF_OUT="${CMAKE_CURRENT_SOURCE_DIR}/results"
  GSH_REPO_DATA="${CMAKE_SOURCE_DIR}/data")
gtest_discover_tests(gsh_kernel_vs_oracle DISCOVERY_TIMEOUT 60 PROPERTIES TIMEOUT 1800)
