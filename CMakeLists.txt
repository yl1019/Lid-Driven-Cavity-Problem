cmake_minimum_required(VERSION 3.1)
project (MainSolver)

add_executable(MainSolver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp ProgramOptions.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(MPI REQUIRED)

target_link_libraries(MainSolver ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${Boost_LIBRARIES})
