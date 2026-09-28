cmake_minimum_required(VERSION 3.10)

project(BaMMmotif VERSION 2.0.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindASan)   # adds the optional build type -DCMAKE_BUILD_TYPE=ASan

# Default to an optimised build. Without a build type CMake compiles with no
# optimisation at all, which makes BaMM!motif about 2x slower.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type given, using Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type: Release, RelWithDebInfo, Debug or ASan" FORCE)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BAMM_NATIVE "Optimise for the CPU of the build machine (-march=native); binaries may not run on other machines" OFF)
option(BAMM_WERROR "Treat compiler warnings as errors" OFF)

add_compile_options(-Wall)
if(BAMM_WERROR)
    add_compile_options(-Werror)
endif()
if(BAMM_NATIVE)
    add_compile_options(-march=native)
endif()

# Boost (header-only parts: special functions and distributions)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# OpenMP is optional; without it all tools run single-threaded
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
    add_compile_definitions(OPENMP)
else()
    message(WARNING "OpenMP not found: BaMM!motif will run single-threaded")
endif()

add_subdirectory(src)

# R and Python helper scripts for plotting and format conversion
install(DIRECTORY R/ DESTINATION bin
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                         GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        FILES_MATCHING PATTERN "*.R")
install(DIRECTORY py/ DESTINATION share/bamm/py
        FILES_MATCHING PATTERN "*.py")

# `ctest` runs all tools on the bundled example and compares the results
# with tests/reference
enable_testing()
add_test(NAME examples
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_examples.sh ${CMAKE_BINARY_DIR}/bin)
