cmake_minimum_required(VERSION 3.11.4)

project(rimerge)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

################################################################################
# Version number
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/.."
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/.."
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hpp")
configure_file(
        ${CMAKE_SOURCE_DIR}/include/rimerge/version.hpp.in
        ${CMAKE_BINARY_DIR}/generated/rimerge/version.hpp)

include_directories(${CMAKE_BINARY_DIR}/generated)

################################################################################
# Compiler and linker flags

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS} ${rimerge-flags}")

################################################################################
# Thirdparty libraries
include(ExternalProject)

# Google Test
ExternalProject_Add(gtest
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "release-1.10.0"
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(gtest INSTALL_DIR)
set(gtest-include ${INSTALL_DIR}/src/gtest/googletest/include)
set(gtest-lib ${INSTALL_DIR}/src/gtest-build/lib/libgtest.a ${INSTALL_DIR}/src/gtest-build/lib/libgtest_main.a)

# SDSL
ExternalProject_Add(sdsl
        GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
        GIT_TAG "v2.1.1"
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl INSTALL_DIR)
set(sdsl-include ${INSTALL_DIR}/src/sdsl-build/include)
set(divsufsort-include ${INSTALL_DIR}/src/sdsl-build/external/libdivsufsort/include)
set(sdsl-lib ${INSTALL_DIR}/src/sdsl-build/lib/libsdsl.a)
set(divsufsort-lib ${INSTALL_DIR}/src/sdsl-build/external/libdivsufsort/lib/libdivsufsort.a)

# spdlog
ExternalProject_Add(spdlog
        GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
        GIT_TAG "v1.5.0"
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}" "-DSPDLOG_MASTER_PROJECT=OFF" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(spdlog INSTALL_DIR)
set(spdlog-include ${INSTALL_DIR}/src/spdlog/include)
set(spdlog-lib ${INSTALL_DIR}/src/spdlog-build/libspdlog.a)

# mio
ExternalProject_Add(mio
        GIT_REPOSITORY "https://github.com/mandreyel/mio.git"
        GIT_TAG ""
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}" "-DBUILD_TESTING=OFF" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(mio INSTALL_DIR)
set(mio-include ${INSTALL_DIR}/src/mio/single_include)

# CLI11
ExternalProject_Add(CLI11
        GIT_REPOSITORY "https://github.com/CLIUtils/CLI11.git"
        GIT_TAG "v1.9.0"
        CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}" "-DCLI11_BUILD_TESTS=OFF" "-DCLI11_BUILD_EXAMPLES=OFF" "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCLI11_BUILD_DOCS=OFF"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(CLI11 INSTALL_DIR)
set(CLI11-include ${INSTALL_DIR}/src/CLI11/include)

# OpenMP, Fixme
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # assumes clang build
        # we can't reliably detect when we're using clang, so for the time being we assume
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(rimerge-includes
        "${CMAKE_SOURCE_DIR}/include"
        "${sdsl-include}"
        "${divsufsort-include}"
        "${spdlog-include}"
        "${mio-include}"
        "${CLI11-include}"
        "${OpenMP_CXX_INCLUDE_DIRS}")
set(rimerge-libs
        "${sdsl-lib}"
        "${divsufsort-lib}"
        "${spdlog-lib}"
        "${OpenMP_CXX_LIBRARY}")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "GNU compiler detected, adding stdc++fs")
    list(APPEND rimerge-libs "stdc++fs")
endif()

################################################################################
# Directory structure

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)


################################################################################
# RImergelib

file (GLOB  ComponentFiles  *  ${PROJECT_SOURCE_DIR}/*.cpp)
add_library(rimerge STATIC ${ComponentFiles})
add_dependencies(rimerge sdsl spdlog mio)
target_include_directories(rimerge PRIVATE "${rimerge-includes}")
target_link_libraries(rimerge ${rimerge_libs})

target_compile_features(rimerge PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Tools
add_subdirectory(tools)

################################################################################
# Tests
link_directories(/usr/local/lib)

add_executable(unit_tests tests/unit_tests.cpp)
target_include_directories(unit_tests PRIVATE "${gtest-include}" "${rimerge-includes}")
add_dependencies(unit_tests rimerge gtest)
target_link_libraries(unit_tests rimerge ${rimerge-libs} ${gtest-lib})

add_executable(bwtmerge_tests tests/test_bwtmerge.cpp)
target_include_directories(bwtmerge_tests PRIVATE "${gtest-include}" "${rimerge-includes}")
add_dependencies(bwtmerge_tests rimerge gtest)
target_link_libraries(bwtmerge_tests rimerge ${rimerge-libs} ${gtest-lib})

################################################################################
# Install


