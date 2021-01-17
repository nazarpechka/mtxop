cmake_minimum_required(VERSION 3.14)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif ()

# Project properties
project(
        "Matrix Operations"
        VERSION 0.1.0
        LANGUAGES CXX)

set(LIBRARY_NAME mtxop)
set(EXECUTABLE_NAME mtxop-test)

# Options
option(BUILD_SHARED_LIBS "Build shared instead of static libraries." OFF)
option(ENABLE_LTO "Enable link time optimization" ON)
option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE
            Release
            CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(
            CACHE CMAKE_BUILD_TYPE
            PROPERTY STRINGS
            "Debug"
            "Release"
            "MinSizeRel"
            "RelWithDebInfo")
endif ()

# Cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(LTO)
include(Warnings)
include(Cache)

# Generate folders for IDE targets (e.g., VisualStudio solutions)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(IDE_FOLDER "")

# Check for LTO support
find_lto(CXX)

# Project warnings
add_library(project_warnings INTERFACE)

set_project_warnings(project_warnings)

# Project options
add_library(project_options INTERFACE)

target_compile_features(project_options INTERFACE cxx_std_17)
# set_target_properties(project_options PROPERTIES 
# 	CXX_STANDARD 17
# 	CXX_STANDARD_REQUIRED YES
# 	CXX_EXTENSIONS NO)

if (MSVC)
    target_compile_options(project_options INTERFACE $<$<CONFIG:RELEASE>:/O2>)
else ()
    target_compile_options(project_options INTERFACE $<$<CONFIG:DEBUG>:-g> $<$<CONFIG:RELEASE>: -funroll-loops -march=native>)
endif ()


# Build library
add_library(${LIBRARY_NAME} src/matrix_ref.cpp
        src/matrix.cpp
        )
target_link_libraries(${LIBRARY_NAME} PRIVATE project_options project_warnings)
target_include_directories(${LIBRARY_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include/${LIBRARY_NAME})


# Build executable
add_executable(${EXECUTABLE_NAME} src/test.cpp)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE project_options project_warnings ${LIBRARY_NAME})
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/external)

target_enable_lto(${EXECUTABLE_NAME} optimized)





