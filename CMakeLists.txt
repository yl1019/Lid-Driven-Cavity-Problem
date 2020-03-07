cmake_minimum_required(VERSION 3.1)
project (MainSolver)

add_executable(MainSolver LidDrivenCavitySolver.cpp LidDrivenCavity.cpp
	ProgramOptions.cpp PoissonSolver.cpp MpiConfig.cpp)

find_package(BLAS REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

target_link_libraries(MainSolver ${BLAS_LIBRARIES} ${Boost_LIBRARIES} ${MPI_LIBRARIES})
