cmake_minimum_required(VERSION 3.12)

# Project Name
project(satutils)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if(TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
endif()

add_library(satutils
  include/satutils/constants.hpp
  src/constants.cpp
  src/binary_ops.cpp
  src/code_gen.cpp
)

target_include_directories(satutils
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(satutils PUBLIC Eigen3::Eigen)
