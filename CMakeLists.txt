cmake_minimum_required(VERSION 3.6)
project(satutils)

string(ASCII 27 Esc)
set(Reset       "${Esc}[m")
set(Bold        "${Esc}[1m")
set(Red         "${Esc}[31m")
set(Green       "${Esc}[32m")
set(Yellow      "${Esc}[33m")
set(Blue        "${Esc}[34m")
set(Magenta     "${Esc}[35m")
set(Cyan        "${Esc}[36m")
set(White       "${Esc}[37m")
set(BoldRed     "${Esc}[1;31m")
set(BoldGreen   "${Esc}[1;32m")
set(BoldYellow  "${Esc}[1;33m")
set(BoldBlue    "${Esc}[1;34m")
set(BoldMagenta "${Esc}[1;35m")
set(BoldCyan    "${Esc}[1;36m")
set(BoldWhite   "${Esc}[1;37m")

# --- Version Control ---
set(SATUTILS_MAJOR_VERSION 0)
set(SATUTILS_MINOR_VERSION 0)
set(SATUTILS_PATCH_VERSION 1)
set(PROJECT_VERSION ${SATUTILS_MAJOR_VERSION}.${SATUTILS_MINOR_VERSION}.${SATUTILS_PATCH_VERSION})
message(STATUS "${BoldBlue}${PROJECT_NAME} version ${PROJECT_VERSION}${Reset}")

# --- Compiler Flags ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wunused -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast -ffast-math")

# --- Packages ---
find_package(Eigen3 REQUIRED)
find_package(navtools REQUIRED)

# --- Add Library Soruce Files ---
set(SATUTILS_HDRS
    include/satutils/atmosphere.hpp
    include/satutils/clock.hpp
    include/satutils/code-gen.hpp
    include/satutils/constants.hpp
    include/satutils/ephemeris.hpp
    include/satutils/parsers.hpp
)

# --- Create the C++ Library ---
add_library(${PROJECT_NAME} INTERFACE ${SATUTILS_HDRS})
target_include_directories(
    ${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${Eigen3_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${Eigen3_INCLUDE_DIRS}>
)
target_link_libraries(
    ${PROJECT_NAME} INTERFACE
    Eigen3::Eigen
    navtools
)

# --- Build Executables ---
if (NOT DEFINED INSTALL_SATUTILS_EXAMPLES OR NOT INSTALL_SATUTILS_EXAMPLES)
    message(STATUS "${BoldCyan}INSTALL_SATUTILS_EXAMPLES is false${Reset}")
else()
    message(STATUS "${BoldCyan}INSTALL_SATUTILS_EXAMPLES is true${Reset}")

    add_executable(test_parser examples/test_parser.cpp)
    target_link_libraries(test_parser PUBLIC ${PROJECT_NAME})
endif()

# --- Make Library 'Findable' for other CMake Packages ---
include(CMakePackageConfigHelpers)

install(
    DIRECTORY include/ 
    DESTINATION include
)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)

export(
    TARGETS ${PROJECT_NAME} 
    FILE ${PROJECT_NAME}Targets.cmake
)

install(
    EXPORT ${PROJECT_NAME}Targets 
    FILE ${PROJECT_NAME}Targets.cmake 
    DESTINATION lib/cmake/${PROJECT_NAME}
)

configure_package_config_file(
    "${PROJECT_NAME}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    INSTALL_DESTINATION lib/cmake/${PROJECT_NAME}
    PATH_VARS PROJECT_NAME
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME}
)
