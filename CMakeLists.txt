cmake_minimum_required(VERSION 3.20)
project(MeshMS LANGUAGES CXX)

# --- C++20 implementation of the MolSurfComp SES algorithm -------------------
# See PORTING_CONTRACT.md (indexing & data-layout conventions) and
# OPTIMIZATION.md (the 28-core, output-preserving parallel design).
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Faithful-first: zero third-party C++ dependencies by default.
option(MESHMS_WITH_CGAL "Build the optional CGAL gate cross-checker (GPL)" OFF)
option(MESHMS_SANITIZE  "Enable ASan/UBSan for the fragile mesher" OFF)
option(MESHMS_NATIVE    "Build with -march=native -O3" OFF)
option(MESHMS_TBB       "Build with Intel oneTBB (parallel stages)" ON)

# Build the CLI / golden+unit tests only when MeshMS is the top-level project
# (a find_package / add_subdirectory consumer should NOT build them). CMake 3.20
# lacks PROJECT_IS_TOP_LEVEL, so detect top-level by source-dir comparison.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(MESHMS_IS_TOP_LEVEL ON)
else()
  set(MESHMS_IS_TOP_LEVEL OFF)
endif()
option(MESHMS_BUILD_TESTS "Build the C++ port unit/golden tests" ${MESHMS_IS_TOP_LEVEL})
option(MESHMS_BUILD_CLI   "Build the MeshMS CLI" ${MESHMS_IS_TOP_LEVEL})

# -ffp-contract=off is MANDATORY: g++ -march=native otherwise fuses a*b+c into
# an FMA, which changes the last-bit rounding and breaks the ULP-sensitive SAS
# golden tests (bit-for-bit 1e-9 topology gate). Keep it always, never optional.
add_compile_options(-Wall -Wextra -ffp-contract=off)
if(MESHMS_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()
if(MESHMS_NATIVE)
  add_compile_options(-march=native -O3)
endif()

# Paths to the xyzr molecule set and the golden reference dumps the tests use.
add_compile_definitions(
  MESHMS_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
  MESHMS_REF_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/ref")

# --- core library: every src/*.cpp (CONFIGURE_DEPENDS = new files auto-picked) -
file(GLOB MESHMS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(MeshMS STATIC ${MESHMS_SOURCES})
# GNUInstallDirs defines CMAKE_INSTALL_INCLUDEDIR/LIBDIR; needed by the
# INSTALL_INTERFACE path below and the install() rules.
include(GNUInstallDirs)
# BUILD_INTERFACE for the in-tree headers; INSTALL_INTERFACE for installed ones
# (required so install(EXPORT ...) can record a relocatable include path).
target_include_directories(MeshMS PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# PRIVATE (not PUBLIC): the library is COMPILED with C++20, but the public facade
# (meshms/capi.hpp) is C++17-clean, so we do NOT force C++20 onto consumers. A
# C++20 static lib links into a C++17 program without issue.
target_compile_features(MeshMS PRIVATE cxx_std_20)

# oneTBB drives the parallel per-entity stages (output-preserving: each entity
# writes its own slot, then a serial ascending-index merge). When OFF -- or when
# TBB is not found -- the code builds and runs SERIALLY and identically:
# meshms/parallel.hpp's serial fallback compiles whenever MESHMS_WITH_TBB is
# undefined, so there is no parallel dependency. TBB comes from the cuemol2
# "deplibs" bundle; a consumer points at it with CMAKE_PREFIX_PATH=<deplibs root>
# (or -DTBB_DIR=<root>/tbb-<ver>/lib/cmake/TBB), which resolves find_package(TBB).
set(MESHMS_TBB_PUBLIC FALSE)  # substituted into the package Config
if(MESHMS_TBB)
  find_package(TBB CONFIG)
  if(TBB_FOUND)
    target_link_libraries(MeshMS PUBLIC TBB::tbb)
    target_compile_definitions(MeshMS PUBLIC MESHMS_WITH_TBB)
    set(MESHMS_TBB_PUBLIC TRUE)
  else()
    message(WARNING "MESHMS_TBB=ON but oneTBB not found; building serially")
  endif()
endif()

# --- CLI ---------------------------------------------------------------------
if(MESHMS_BUILD_CLI AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/cli/main.cpp)
  add_executable(meshms_cli ${CMAKE_CURRENT_SOURCE_DIR}/src/cli/main.cpp)
  target_link_libraries(meshms_cli PRIVATE MeshMS)
endif()

# --- tests: one executable per tests/test_*.cpp, registered with CTest --------
if(MESHMS_BUILD_TESTS)
  enable_testing()
  file(GLOB TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp)
  foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE MeshMS)
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()
endif()

# --- install + find_package(MeshMS) ------------------------------------------
# Install the static lib + ALL public headers, export a namespaced target, and
# generate a relocatable package Config so an external project can do
# `find_package(MeshMS REQUIRED)` and link `MeshMS::MeshMS`.
include(CMakePackageConfigHelpers)

install(TARGETS MeshMS
        EXPORT MeshMSTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# include/ holds the meshms/ subdir; installs to <prefix>/include/meshms/*.hpp.
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT MeshMSTargets
        FILE MeshMSTargets.cmake
        NAMESPACE MeshMS::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshMS)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MeshMSConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MeshMSConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshMS)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MeshMSConfigVersion.cmake
  VERSION 0.1.0
  COMPATIBILITY SameMajorVersion)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/MeshMSConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/MeshMSConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshMS)
