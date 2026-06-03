cmake_minimum_required(VERSION 3.26)

cmake_policy(SET CMP0057 NEW) # ``IN_LIST will be interpreted as an operator when the policy is set to NEW``
cmake_policy(SET CMP0074 NEW) # ``find_package()`` uses ``<PackageName>_ROOT`` variables.
cmake_policy(SET CMP0077 NEW) # ``option()`` honors normal variables.
cmake_policy(SET CMP0080 NEW) # ``BundleUtilities`` cannot be included at configure time.

# The ``FindPythonInterp`` and ``FindPythonLibs`` modules are removed.
if(CMAKE_VERSION VERSION_GREATER "3.27.0")
  cmake_policy(SET CMP0148 NEW)
endif()

# resolves symlinks before collapsing ../ components.
if(CMAKE_VERSION VERSION_GREATER "3.28.0")
  cmake_policy(SET CMP0152 NEW)
endif()

if(CMAKE_VERSION VERSION_GREATER "3.29.0")
  cmake_policy(SET CMP0167 NEW) # Suppress FindBoost deprecation warning (retained for EbsdLib compatibility)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmake/Utility.cmake)

include(CMakeDependentOption)


if(${MTRSIM_BUILD_STANDALONE_LIB})
  message(STATUS "MTRSim as Library")
  project(MTRSim
    VERSION 1.0.0
    DESCRIPTION "Microtexture Region Simulator"
    LANGUAGES CXX
  )
  include(LibMTRSim.cmake)
else()
  message(STATUS "MTRSim as SimplnxPlugin")
  include(MTRSimPlugin.cmake)
endif()