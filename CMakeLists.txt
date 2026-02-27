cmake_minimum_required(VERSION 3.26)

cmake_policy(SET CMP0074 NEW) # ``find_package()`` uses ``<PackageName>_ROOT`` variables.
cmake_policy(SET CMP0077 NEW) # ``option()`` honors normal variables.
cmake_policy(SET CMP0080 NEW) # ``BundleUtilities`` cannot be included at configure time.
cmake_policy(SET CMP0167 NEW) # ``The FindBoost module is removed.

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

project(MTRsim
  VERSION 1.0.0
  DESCRIPTION "Microtexture Region Simulator"
  LANGUAGES CXX
)

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
# Dependencies
# ------------------------------------------------------------------------------
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS math)
find_package(CLI11 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(Stb REQUIRED)
find_package(HDF5 1.14 MODULE REQUIRED)
get_target_property(hdf5_dll_path hdf5::hdf5 IMPORTED_LOCATION_RELEASE)
get_filename_component(hdf5_dll_path "${hdf5_dll_path}" DIRECTORY)
get_property(MTRSIM_EXTRA_LIBRARY_DIRS GLOBAL PROPERTY MTRSIM_EXTRA_LIBRARY_DIRS)
set_property(GLOBAL PROPERTY MTRSIM_EXTRA_LIBRARY_DIRS ${MTRSIM_EXTRA_LIBRARY_DIRS} ${hdf5_dll_path})

# if(NOT TARGET H5Support::H5Support)
#   add_subdirectory("${MTRsim_SOURCE_DIR}/../H5Support" "${MTRsim_BINARY_DIR}/H5Support")
# endif()

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


if(NOT TARGET EbsdLib::EbsdLib)
    if(EXISTS "${MTRsim_SOURCE_DIR}/../EbsdLib")
        set(EbsdLibProj_SOURCE_DIR "${MTRsim_SOURCE_DIR}/../EbsdLib")
    else()
        message(FATAL_ERROR "EbsdLibProj_SOURCE_DIR was not set. Where is the EbsdLib project directory. Please set the EbsdLibProj_SOURCE_DIR variable to the EbsdLib directory.")
    endif()
    message(STATUS "EbsdLibProj_SOURCE_DIR: ${EbsdLibProj_SOURCE_DIR}")

    set(EbsdLib_ENABLE_HDF5 OFF)
    set(EbsdLib_USE_PARALLEL_ALGORITHMS ${MTRSIM_ENABLE_MULTICORE})
    set(EbsdLib_BUILD_H5SUPPORT OFF)
    set(H5Support_INCLUDE_QT_API OFF)
    set(EbsdLib_BUILD_TESTS OFF)
    add_subdirectory( ${EbsdLibProj_SOURCE_DIR} ${PROJECT_BINARY_DIR}/EbsdLib)
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

if(MTRSIM_BUILD_TESTS)
  add_subdirectory(tests)
endif()
