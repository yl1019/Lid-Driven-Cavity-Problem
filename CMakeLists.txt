cmake_minimum_required(VERSION 2.8.9)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

project (MainSolver)
add_executable(MainSolver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(MainSolver ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
