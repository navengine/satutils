cmake_minimum_required(VERSION 3.5)

# Project Name
project(satutils)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(satutils
  include/satutils/constants.hpp
  src/constants.cpp
  src/binary_ops.cpp
  src/code_gen.cpp
)

target_include_directories(satutils
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

