cmake_minimum_required(VERSION 3.10)

project(RyadalevskayaSolver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#include(CheckCXXCompilerFlags)
#include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTS "Build tests" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

if(Eigen3_FOUND)
	message(STATUS "Eigen found: ${EIGEN3_VERSION}")
	message(STATUS "Eigen include dir: ${EIGEN3_INCLUDE_DIR}")
else()
	message(FATAL_ERROR "Eigen not found!")
endif()

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${EIGEN3_INCLUDE_DIR}
)

add_executable(ryad-solver
	src/core/main.cpp
	src/numerics/solver.cpp
	src/thermo/statSum.cpp
)

target_include_directories(ryad-solver PRIVATE ${EIGEN3_INCLUDE_DIR})


