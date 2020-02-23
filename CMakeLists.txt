cmake_minimum_required(VERSION 2.8.9)
project (MainSolver)
add_executable(MainSolver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(MainSolver ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
