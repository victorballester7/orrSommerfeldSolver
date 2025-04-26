cmake_minimum_required(VERSION 3.10)
project(OrrSommerfeldSolver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen package
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Add executable
add_executable(orr_sommerfeld_solver orr_sommerfeld_solver.cpp)
