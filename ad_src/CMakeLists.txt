# ----------------------------------------------------------------------
# Standalone build for the ad_src/ drag(alpha) drivers.
#
# This is a self-contained CMake project, deliberately NOT wired into
# sparta/cmake or sparta/src/CMakeLists.txt (no add_subdirectory anywhere
# references ad_src). Configuring/building the main SPARTA CMake tree
# (sparta/build) is completely unaffected by this file.
#
# SPARTA itself is built as a side effect, via CMake custom commands that
# shell out to SPARTA's own classic src/Makefile in "mode=lib" (the same
# mechanism `cd src && make mode=lib serial` uses upstream): it produces
# src/libsparta_serial.a and, if DRAG_ENABLE_AD is on, src/libsparta_serial_ad.a
# plus src/STUBS/libmpi_stubs.a. Those live in sparta/src/, alongside (but
# not overwriting) any per-machine .o dirs (Obj_serial, Obj_serial_ad) the
# classic Makefile already uses -- nothing here touches sparta/build.
#
# Usage, from ad_src/:
#   cmake -S . -B build                     # stock targets only
#   cmake -S . -B build -DDRAG_ENABLE_AD=ON  # also build the AD variant
#   cmake --build build
#   ./build/drag_run --input drag.in         # run from a dir with air.species/air.vss
#
# cmake -S . -B build -DHELP=ON   prints this usage and stops, without
# requiring a C++ compiler or attempting to build anything.
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16)
project(sparta_drag_tools NONE)

option(HELP "Print build instructions and exit, without configuring a compiler" OFF)
if(HELP)
  message(NOTICE
    "ad_src build instructions\n"
    "---------------------------------------------------------------\n"
    "SPARTA does NOT need to be pre-built. Configuring/building this\n"
    "project builds SPARTA itself as a side effect, via SPARTA's own\n"
    "classic src/Makefile in 'mode=lib' (static-library mode), scoped\n"
    "to the 'serial' (MPI-stubs, no real MPI) machine target. Output\n"
    "libraries land in sparta/src/ (libsparta_serial.a, libmpi_stubs.a);\n"
    "sparta/build (the separate CMake build of SPARTA, if any) is never\n"
    "touched.\n"
    "\n"
    "Stock build (drag_run, opt_run, bezier_gen):\n"
    "    cmake -S . -B build\n"
    "    cmake --build build\n"
    "\n"
    "AD build (adds drag_run_ad, test_ad_geom, built against a second\n"
    "SPARTA library compiled with -DSPARTA_AD -- the forward-mode dual\n"
    "number sfloat type, see src/sfloat.h). This recompiles all of\n"
    "SPARTA a second time under a separate object dir (Obj_serial_ad),\n"
    "so first configure is slower; stock targets are unaffected:\n"
    "    cmake -S . -B build -DDRAG_ENABLE_AD=ON\n"
    "    cmake --build build\n"
    "\n"
    "Run drivers from a directory containing air.species/air.vss, e.g.\n"
    "bezier. bezier_gen has no SPARTA dependency.\n"
    "---------------------------------------------------------------")
  return()
endif()

enable_language(CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(DRAG_ENABLE_AD "Also build drag_run_ad / test_ad_geom against a -DSPARTA_AD build of SPARTA" OFF)

set(SPARTA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(SPARTA_SRC_DIR ${SPARTA_ROOT}/src)
find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

# ---- STUBS (MPI stand-in for serial builds), shared by both variants --

set(SPARTA_MPI_STUBS_A ${SPARTA_SRC_DIR}/STUBS/libmpi_stubs.a)
add_custom_command(
  OUTPUT ${SPARTA_MPI_STUBS_A}
  COMMAND ${MAKE_EXECUTABLE} -C ${SPARTA_SRC_DIR}/STUBS
  COMMENT "Building SPARTA MPI stubs (src/STUBS)")
add_custom_target(sparta_build_mpi_stubs DEPENDS ${SPARTA_MPI_STUBS_A})
add_library(mpi_stubs STATIC IMPORTED)
set_target_properties(mpi_stubs PROPERTIES IMPORTED_LOCATION ${SPARTA_MPI_STUBS_A})
add_dependencies(mpi_stubs sparta_build_mpi_stubs)

# ---- stock SPARTA library: classic `make mode=lib serial` -------------

set(SPARTA_SERIAL_A ${SPARTA_SRC_DIR}/libsparta_serial.a)
add_custom_command(
  OUTPUT ${SPARTA_SERIAL_A}
  COMMAND ${MAKE_EXECUTABLE} -C ${SPARTA_SRC_DIR} mode=lib serial
  DEPENDS ${SPARTA_MPI_STUBS_A}
  COMMENT "Building SPARTA (serial, static lib)")
add_custom_target(sparta_build_serial DEPENDS ${SPARTA_SERIAL_A})
add_library(sparta_serial STATIC IMPORTED)
set_target_properties(sparta_serial PROPERTIES IMPORTED_LOCATION ${SPARTA_SERIAL_A})
add_dependencies(sparta_serial sparta_build_serial)

# drag_run: links against the SPARTA library through its public C
# interface (src/library.h).
add_executable(drag_run drag_main.cpp drag_objective.cpp bezier_geom.cpp)
target_compile_definitions(drag_run PRIVATE SPARTA_UNORDERED_MAP)
target_include_directories(drag_run PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
target_link_libraries(drag_run PRIVATE sparta_serial mpi_stubs)

# opt_run: same SPARTA dependency, different driver (optimization loop).
add_executable(opt_run opt_main.cpp drag_objective.cpp bezier_geom.cpp)
target_compile_definitions(opt_run PRIVATE SPARTA_UNORDERED_MAP)
target_include_directories(opt_run PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
target_link_libraries(opt_run PRIVATE sparta_serial mpi_stubs)

# bezier_gen: standalone geometry generator, no SPARTA dependency.
add_executable(bezier_gen bezier_gen.cpp bezier_geom.cpp)

# ---- AD SPARTA library: `make mode=lib serial_ad` (-DSPARTA_AD) -------

if(DRAG_ENABLE_AD)
  set(SPARTA_SERIAL_AD_A ${SPARTA_SRC_DIR}/libsparta_serial_ad.a)
  add_custom_command(
    OUTPUT ${SPARTA_SERIAL_AD_A}
    COMMAND ${MAKE_EXECUTABLE} -C ${SPARTA_SRC_DIR} mode=lib serial_ad
    DEPENDS ${SPARTA_MPI_STUBS_A}
    COMMENT "Building SPARTA (serial, -DSPARTA_AD, static lib)")
  add_custom_target(sparta_build_serial_ad DEPENDS ${SPARTA_SERIAL_AD_A})
  add_library(sparta_serial_ad STATIC IMPORTED)
  set_target_properties(sparta_serial_ad PROPERTIES IMPORTED_LOCATION ${SPARTA_SERIAL_AD_A})
  add_dependencies(sparta_serial_ad sparta_build_serial_ad)

  # drag_run_ad: same sources as drag_run, compiled with -DSPARTA_AD and
  # linked against the AD SPARTA library instead. bezier_geom.cpp has no
  # AD-conditional code, so it's built plain and shared with the stock
  # targets above (matches the old Makefile's behavior).
  add_executable(drag_run_ad drag_main.cpp drag_objective.cpp bezier_geom.cpp)
  target_compile_definitions(drag_run_ad PRIVATE SPARTA_AD SPARTA_UNORDERED_MAP)
  target_include_directories(drag_run_ad PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
  target_link_libraries(drag_run_ad PRIVATE sparta_serial_ad mpi_stubs)

  add_executable(test_ad_geom test_ad_geom.cpp bezier_geom.cpp)
  target_compile_definitions(test_ad_geom PRIVATE SPARTA_AD SPARTA_UNORDERED_MAP)
  target_include_directories(test_ad_geom PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
  target_link_libraries(test_ad_geom PRIVATE sparta_serial_ad mpi_stubs)
endif()
