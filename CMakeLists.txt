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
option(MESHMS_LTO       "Enable interprocedural optimization (LTO)" OFF)
# Target ISA baseline for redistributable builds; see the MESHMS_FP block below.
set(MESHMS_ARCH "" CACHE STRING "ISA baseline: 'avx2' (portable alias) or a raw -march=/arch: value")

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
option(MESHMS_BUILD_TOOLS "Build the bench/converter tools" ${MESHMS_IS_TOP_LEVEL})

# --- floating-point policy (MESHMS_FP) ---------------------------------------
# Selects the FP semantics the library is compiled with.
#
#   strict (default) -- the bit-exact build the golden suite gates on. NO FMA
#     contraction: a fused a*b+c changes the last-bit rounding and flips the
#     ULP-sensitive SAS boundary tests, which changes the triple-point set and
#     with it the integer topology (see OPTIMIZATION.md, PORTING_CONTRACT.md).
#     pysq() stays pow(x, 2.0). Use this for development, bug reproduction, and
#     anything compared against the golden fixtures.
#
#   fast -- deploy builds (cuemol2/cuemol3). Trades the bit-for-bit golden gate
#     for speed: FMA contraction, no math errno, no trapping math, no signed
#     zeros, reciprocal division, and pysq() == x*x. Verified by the equivalence
#     gate (tests/test_fp_gate.cpp) instead of the golden gate.
#
# Deliberately NOT part of fast, and rejected outright below:
#   -ffast-math / -Ofast / -ffinite-math-only
#     -ffinite-math-only folds isfinite()/isnan() to a constant, which disables
#     the NaN detection the deploy gate rests on.
#   -ffast-math / -Ofast / -funsafe-math-optimizations
#     all three match GCC's crtfastmath.o link spec, which sets FTZ/DAZ for the
#     WHOLE host process -- unacceptable for a library linked into cuemol2. The
#     individual flags fast does use match no such spec.
#   -fassociative-math (implied by -funsafe-math-optimizations)
#     the parallel design deliberately has no float reductions, so reassociation
#     buys almost nothing, while it WOULD reassociate the sequential sums in
#     mesh_area/signed_volume -- the very quantities the gate measures.
set(MESHMS_FP "strict" CACHE STRING "FP policy: strict (bit-exact golden gate) | fast (deploy)")
set_property(CACHE MESHMS_FP PROPERTY STRINGS strict fast)
if(NOT MESHMS_FP MATCHES "^(strict|fast)$")
  message(FATAL_ERROR "MESHMS_FP must be 'strict' or 'fast' (got '${MESHMS_FP}')")
endif()

# Reject the unsafe FP flags in BOTH policies, however they were injected. The
# header guard in meshms/vec3.hpp catches -ffast-math even when it arrives from
# outside CMake; this check reports it at configure time with a usable message.
option(MESHMS_ALLOW_UNSAFE_FP_FLAGS "Skip the -ffast-math/-Ofast rejection (at your own risk)" OFF)
if(NOT MESHMS_ALLOW_UNSAFE_FP_FLAGS)
  set(_meshms_fp_scan "${CMAKE_CXX_FLAGS}")
  foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    string(APPEND _meshms_fp_scan " ${CMAKE_CXX_FLAGS_${_cfg}}")
  endforeach()
  foreach(_bad "-ffast-math" "-Ofast" "-ffinite-math-only" "-funsafe-math-optimizations"
               "-fapprox-func" "/fp:fast")
    string(FIND "${_meshms_fp_scan}" "${_bad}" _hit)
    if(NOT _hit EQUAL -1)
      message(FATAL_ERROR
              "${_bad} is not supported by MeshMS: it either disables isfinite() "
              "or sets FTZ/DAZ for the whole host process. Use -DMESHMS_FP=fast "
              "for a fast deploy build, or -DMESHMS_ALLOW_UNSAFE_FP_FLAGS=ON to "
              "override at your own risk.")
    endif()
  endforeach()
endif()

# MESHMS_ARCH takes the raw value for the compiler in use (-march=<v> on
# GNU/Clang, /arch:<v> on MSVC), plus one portable alias:
#
#   avx2 -- the x86-64 deploy baseline. Expands to -march=x86-64-v3 (AVX2 + FMA +
#           BMI) or /arch:AVX2, and is IGNORED on non-x86 targets, so one build
#           script can pass it on every platform. FMA is the point: with
#           MESHMS_FP=fast, neither /fp:contract nor -ffp-contract=fast can emit
#           an FMA unless the target actually has the instruction, and the
#           default x86-64 baseline does not.
set(_meshms_arch "${MESHMS_ARCH}")
if(_meshms_arch STREQUAL "avx2")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^([xX]86_64|AMD64|amd64)$")
    if(MSVC)
      set(_meshms_arch "AVX2")
    else()
      set(_meshms_arch "x86-64-v3")
    endif()
  else()
    message(STATUS "MeshMS: MESHMS_ARCH=avx2 ignored on ${CMAKE_SYSTEM_PROCESSOR}")
    set(_meshms_arch "")
  endif()
endif()

if(MSVC)
  add_compile_options(/W4 /utf-8 /fp:precise)
  if(MSVC_VERSION GREATER_EQUAL 1930)
    if(MESHMS_FP STREQUAL "strict")
      add_compile_options(/fp:contract:off)
    else()
      # fast stops at contraction on MSVC: /fp:fast has no finer granularity and
      # would also relax the NaN/Inf assumptions the gate depends on. On MSVC the
      # win comes from pysq() == x*x instead -- /fp:precise never folds pow(x, 2.0).
      # Note FMA only materialises with /arch:AVX2 or later (see MESHMS_ARCH).
      add_compile_options(/fp:contract)
    endif()
  endif()
  if(_meshms_arch)
    add_compile_options(/arch:${_meshms_arch})
  endif()
  # MESHMS_SANITIZE / MESHMS_NATIVE stay GNU/Clang-only.
else()
  add_compile_options(-Wall -Wextra)
  if(MESHMS_FP STREQUAL "strict")
    add_compile_options(-ffp-contract=off)
  else()
    add_compile_options(-ffp-contract=fast -fno-math-errno -fno-trapping-math
                        -fno-signed-zeros -freciprocal-math)
  endif()
  if(MESHMS_SANITIZE)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
  endif()
  if(MESHMS_NATIVE)
    add_compile_options(-march=native -O3)
  endif()
  if(_meshms_arch)
    # Explicit target baseline for redistributable builds: -march=native bakes in
    # the build machine's ISA and SIGILLs elsewhere.
    add_compile_options(-march=${_meshms_arch})
  endif()
endif()

# LTO is opt-in and global (library + CLI + tools + tests, so a benchmark build
# measures the whole pipeline). A static library built with LTO carries bitcode,
# so a consumer linking it must use a matching toolchain.
if(MESHMS_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT MESHMS_IPO_OK OUTPUT MESHMS_IPO_MSG)
  if(MESHMS_IPO_OK)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "MESHMS_LTO=ON but IPO is unavailable: ${MESHMS_IPO_MSG}")
  endif()
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

# A non-strict FP policy drops the bit-for-bit golden gate, so the tests switch
# to the equivalence gate (test_fp_gate) and version() reports the policy. Not
# defined for strict, so a default build is unchanged down to the version string.
if(NOT MESHMS_FP STREQUAL "strict")
  target_compile_definitions(MeshMS PUBLIC MESHMS_FP_FAST MESHMS_FP_MODE="${MESHMS_FP}")
endif()

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

# --- tools: per-stage benchmark + xyz/PDB -> xyzr converter -------------------
# Not installed; for local measurement work (see OPTIMIZATION.md).
if(MESHMS_BUILD_TOOLS)
  add_executable(meshms_bench ${CMAKE_CURRENT_SOURCE_DIR}/tools/bench.cpp)
  target_link_libraries(meshms_bench PRIVATE MeshMS)
  add_executable(meshms_mol2xyzr ${CMAKE_CURRENT_SOURCE_DIR}/tools/mol2xyzr.cpp)
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
    # Golden bit-regression tests return 77 under MESHMS_FP=fast (see
    # MESHMS_SKIP_IF_FAST in tests/test_util.hpp) so ctest reports them as
    # "Skipped" rather than silently passing. Harmless in a strict build.
    set_tests_properties(${test_name} PROPERTIES SKIP_RETURN_CODE 77)
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
