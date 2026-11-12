find_package(GTest CONFIG REQUIRED)
add_executable(gsh_tests
  test_legendre.cpp
  test_model.cpp
  test_shs.cpp
  test_functionals.cpp
  test_grid.cpp
  test_protocol.cpp
  test_tides.cpp
  test_coordinates.cpp
  test_errors.cpp
  test_error_report.cpp
  test_sanity_identities.cpp
)
target_link_libraries(gsh_tests PRIVATE gsh GTest::gtest GTest::gtest_main)
target_include_directories(gsh_tests PRIVATE ${CMAKE_SOURCE_DIR}/src/cpp/gsh_legendre/kernels)
# GSH_REPO_DATA points at the repository data folder, which the public tree does not have. Tests that need
# external model or grid files (EGM2008, XGM2019e_2159, Colorado grid) skip themselves when the file is absent.
target_compile_definitions(gsh_tests PRIVATE
  GSH_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
  GSH_REPO_DATA="${CMAKE_SOURCE_DIR}/data")
include(GoogleTest)
gtest_discover_tests(gsh_tests DISCOVERY_TIMEOUT 60)
