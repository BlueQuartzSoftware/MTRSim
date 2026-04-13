cmake_minimum_required(VERSION 3.26)

cmake_policy(SET CMP0074 NEW) # ``find_package()`` uses ``<PackageName>_ROOT`` variables.
cmake_policy(SET CMP0077 NEW) # ``option()`` honors normal variables.
cmake_policy(SET CMP0080 NEW) # ``BundleUtilities`` cannot be included at configure time.
cmake_policy(SET CMP0167 NEW) # Suppress FindBoost deprecation warning (retained for EbsdLib compatibility)

# The ``FindPythonInterp`` and ``FindPythonLibs`` modules are removed.
# if(CMAKE_VERSION VERSION_GREATER "3.27.0")
#   cmake_policy(SET CMP0148 NEW)
# endif()

# resolves symlinks before collapsing ../ components.
if(CMAKE_VERSION VERSION_GREATER "3.28.0")
  cmake_policy(SET CMP0152 NEW)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmake/Utility.cmake)


include(CMakeDependentOption)

project(MTRSim
  VERSION 1.0.0
  DESCRIPTION "Microtexture Region Simulator"
  LANGUAGES CXX
)

set(MTRSim_BIN_DIR ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${MTRSim_BIN_DIR} CACHE PATH "Single Directory for all Libraries")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MTRSim_BIN_DIR} CACHE PATH "Single Directory for all Executables.")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${MTRSim_BIN_DIR} CACHE PATH "Single Directory for all static libraries.")

list(APPEND CMAKE_MODULE_PATH ${MTRSim_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# ------------------------------------------------------------------------------
# is building unit tests enabled
# ------------------------------------------------------------------------------
option(MTRSIM_BUILD_TESTS "Enable buildingtests" ON)
enable_vcpkg_manifest_feature(TEST_VAR MTRSIM_BUILD_TESTS FEATURE "tests")

# ------------------------------------------------------------------------------
# are multithreading algorithms enabled
# ------------------------------------------------------------------------------
option(MTRSIM_ENABLE_MULTICORE "Enable multicore support" ON)
enable_vcpkg_manifest_feature(TEST_VAR MTRSIM_ENABLE_MULTICORE FEATURE "parallel")

# ------------------------------------------------------------------------------
# Python bindings (optional)
# ------------------------------------------------------------------------------
option(MTRSIM_BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" ON)
enable_vcpkg_manifest_feature(TEST_VAR MTRSIM_BUILD_PYTHON_BINDINGS FEATURE "python")

if(MTRSIM_BUILD_PYTHON_BINDINGS)
  find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
endif()

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
find_package(Eigen3 CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(Stb REQUIRED)

# -----------------------------------------------------------------------
# Find HDF5 and get the path to the DLL libraries and put that into a
# global property for later install, debugging and packaging
# -----------------------------------------------------------------------
find_package(HDF5 1.14 MODULE REQUIRED)
get_target_property(hdf5_dll_path hdf5::hdf5 IMPORTED_LOCATION_RELEASE)
get_filename_component(hdf5_dll_path "${hdf5_dll_path}" DIRECTORY)
get_property(MTRSim_EXTRA_LIBRARY_DIRS GLOBAL PROPERTY MTRSim_EXTRA_LIBRARY_DIRS)
set_property(GLOBAL PROPERTY MTRSim_EXTRA_LIBRARY_DIRS ${MTRSim_EXTRA_LIBRARY_DIRS} ${hdf5_dll_path})

# -----------------------------------------------------------------------
# Find oneTBB and get the path to the DLL libraries and put that into a
# global property for later install, debugging and packaging
# -----------------------------------------------------------------------
if(MTRSIM_ENABLE_MULTICORE)
  find_package(TBB CONFIG REQUIRED)
  get_target_property(tbb_dll_path TBB::tbb IMPORTED_LOCATION_RELEASE)
  get_filename_component(tbb_dll_path "${tbb_dll_path}" DIRECTORY)
  get_property(MTRSIM_EXTRA_LIBRARY_DIRS GLOBAL PROPERTY MTRSIM_EXTRA_LIBRARY_DIRS)
  set_property(GLOBAL PROPERTY MTRSIM_EXTRA_LIBRARY_DIRS ${MTRSIM_EXTRA_LIBRARY_DIRS} ${tbb_dll_path})
endif()


# ------------------------------------------------------------------------------
# Is the OrientationAnalysis Plugin enabled [DEFAULT=ON]
# ------------------------------------------------------------------------------
option(MTRSIM_ENABLE_OrientationAnalysis "Enable the OrientationAnalysis Plugin" ON)
if(NOT MTRSIM_USE_LOCAL_EBSD_LIB)
  enable_vcpkg_manifest_feature(TEST_VAR MTRSIM_ENABLE_OrientationAnalysis FEATURE "ebsd")
endif()

# ------------------------------------------------------------------------------
# Required EbsdLib and H5Support
# ------------------------------------------------------------------------------
if(MTRSIM_USE_LOCAL_EBSD_LIB)
  include("${MTRSim_SOURCE_DIR}/cmake/UseEbsdLib.cmake")
else()
  find_package(H5Support REQUIRED)
  find_package(EbsdLib REQUIRED)
endif()


if(MTRSIM_BUILD_TESTS)
  find_package(Catch2 CONFIG REQUIRED)
  include(CTest)
  include(Catch)
endif()

# ------------------------------------------------------------------------------
# Sub-projects
# ------------------------------------------------------------------------------
add_subdirectory(src/libmtrsim)
add_subdirectory(app)
add_subdirectory(tools)

if(MTRSIM_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(MTRSIM_BUILD_PYTHON_BINDINGS)
  add_subdirectory(wrapping/python)
endif()
