add_executable(ex_puntual ex_puntual.cpp)
target_link_libraries(ex_puntual PRIVATE gsh)

add_executable(ex_legendre ex_legendre.cpp)
target_link_libraries(ex_legendre PRIVATE gsh)

add_executable(ex_functionals_table ex_functionals_table.cpp)
target_link_libraries(ex_functionals_table PRIVATE gsh)

add_executable(ex_grid ex_grid.cpp)
target_link_libraries(ex_grid PRIVATE gsh)

add_executable(ex_errors ex_errors.cpp)
target_link_libraries(ex_errors PRIVATE gsh)

add_executable(ex_ellipsoid_tides ex_ellipsoid_tides.cpp)
target_link_libraries(ex_ellipsoid_tides PRIVATE gsh)

add_executable(ex_model_inspect ex_model_inspect.cpp)
target_link_libraries(ex_model_inspect PRIVATE gsh)
