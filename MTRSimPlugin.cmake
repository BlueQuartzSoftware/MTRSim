# ------------------------------------------------------------------------------
# Required EbsdLib and H5Support
# ------------------------------------------------------------------------------
if(SIMPLNX_USE_LOCAL_EBSD_LIB)
  include("${simplnx_SOURCE_DIR}/cmake/UseEbsdLib.cmake")
else()
  find_package(H5Support REQUIRED)
  find_package(EbsdLib REQUIRED)
endif()

# ------------------------------------------------------------------------------
# EbsdLib needs install rules for creating packages
# ------------------------------------------------------------------------------
get_property(GENERATOR_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

#if(NOT GENERATOR_IS_MULTI_CONFIG)
  get_property(SIMPLNX_EXTRA_LIBRARY_DIRS GLOBAL PROPERTY SIMPLNX_EXTRA_LIBRARY_DIRS)
  set_property(GLOBAL PROPERTY SIMPLNX_EXTRA_LIBRARY_DIRS ${SIMPLNX_EXTRA_LIBRARY_DIRS} ${EbsdLib_LIB_DIRS})
#endif()

# ------------------------------------------------------------------------------
# Include the standard simplnx plugin cmake support codes
# ------------------------------------------------------------------------------
include("${simplnx_SOURCE_DIR}/cmake/Plugin.cmake")

set(PLUGIN_NAME "MTRSim")
set(${PLUGIN_NAME}_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

# ------------------------------------------------------------------------------
# These are all the filters in the plugin. All filters should be kept in the
# MTRSim/src/MTRSim/Filters/ directory.
set(FilterList
    ComputeODFFilter
    ReadMTRSimODFFilter
    WriteMTRSimODFFilter
)
  
set(ActionList
)

# ------------------------------------------------------------------------------
# Algorithms
# This should be integrated with the `create_simplnx_plugin` function call
# ------------------------------------------------------------------------------
set(AlgorithmList
    ComputeODF
    ReadMTRSimODF
    WriteMTRSimODF
)

create_simplnx_plugin(NAME ${PLUGIN_NAME}
  FILTER_LIST ${FilterList}
  ACTION_LIST ${ActionList}
  ALGORITHM_LIST ${AlgorithmList}  
  DESCRIPTION "${PLUGIN_NAME} Plugin"
  VERSION "1.0.0"
  DOC_CHECK
  ADD_TO_GLOBAL_LIST
)

#------------------------------------------------------------------------------
# Add in the unit tests for the filters included in this plugin
if(SIMPLNX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(${${PLUGIN_NAME}_SOURCE_DIR}/test)
endif()



# ------------------------------------------------------------------------------
# Add any filter parameters to the `simplnx` target
# This should be integrated with the `create_simplnx_plugin` function call
#
# The `simplnx` target will need to know about any dependent libraries and any
# include directories that are needed for the parameter.
# ------------------------------------------------------------------------------
set(simplnx_injected_parameters

)

foreach(source ${simplnx_injected_parameters})
  target_sources(simplnx PRIVATE
    ${${PLUGIN_NAME}_SOURCE_DIR}/src/${PLUGIN_NAME}/Parameters/${source}.h
    ${${PLUGIN_NAME}_SOURCE_DIR}/src/${PLUGIN_NAME}/Parameters/${source}.cpp)

  source_group(TREE "${${PLUGIN_NAME}_SOURCE_DIR}/src/${PLUGIN_NAME}" PREFIX ${PLUGIN_NAME}
    FILES ${${PLUGIN_NAME}_SOURCE_DIR}/src/${PLUGIN_NAME}/Parameters/${source}.h
    ${${PLUGIN_NAME}_SOURCE_DIR}/src/${PLUGIN_NAME}/Parameters/${source}.cpp)
endforeach()

target_include_directories(simplnx PUBLIC )

# ------------------------------------------------------------------------------
# Add a link library to simplnx because it now depends on EbsdLib. Simplnx will
# be able to find the include dirs from the CMake EbsdLib target
# ------------------------------------------------------------------------------
#target_link_libraries(simplnx PUBLIC EbsdLib::EbsdLib)


# ------------------------------------------------------------------------------
# Extra Sources
# ------------------------------------------------------------------------------
set(PLUGIN_EXTRA_SOURCES
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/AssignmentRule.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/GPGenerator.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/IPFMapper.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/MTRDataLoader.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFBuilder.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFCalculator.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFFileIO.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFSampler.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/PGRFSimulation.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/PoleFigure.cpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/QSimVN.cpp
)
set(PLUGIN_EXTRA_HEADERS
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/AssignmentRule.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/GPGenerator.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/IPFMapper.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/MTRDataLoader.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFBuilder.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFCalculator.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFFileIO.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/ODFSampler.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/PGRFSimulation.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/PoleFigure.hpp
  ${${PLUGIN_NAME}_SOURCE_DIR}/src/LibMTRSim/QSimVN.hpp
)

target_sources(${PLUGIN_NAME} 
  PRIVATE
    ${PLUGIN_EXTRA_SOURCES}
    ${PLUGIN_EXTRA_HEADERS}
)

#------------------------------------------------------------------------------
# If there are additional library packages that need to be found, this is where
# to do that using the usual find_package(....) call
# find_package(....)


# ------------------------------------------------------------------------------
# If there are additional libraries that this plugin needs to link against you
# can use the target_link_libraries() cmake call
# ------------------------------------------------------------------------------
target_link_libraries(${PLUGIN_NAME} PUBLIC EbsdLib::EbsdLib)

# ------------------------------------------------------------------------------
# If there are additional include directories that are needed for this plugin
# you can use the target_include_directories(.....) cmake call
target_include_directories(${PLUGIN_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/MTRSim>

)

# -----------------------------------------------------------------------
# Install example pipelines
# -----------------------------------------------------------------------
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(is_multi_config)
  set(CX_CONFIG_DIR "$<CONFIG>")
else()
  set(CX_CONFIG_DIR ".")
endif()

if(EXISTS "${${PLUGIN_NAME}_SOURCE_DIR}/pipelines")
  set(PIPELINE_DEST_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CX_CONFIG_DIR}/pipelines/${PLUGIN_NAME}")
  add_custom_target(Copy_${PLUGIN_NAME}_Pipeline_Folder ALL
                    COMMAND ${CMAKE_COMMAND} -E copy_directory "${${PLUGIN_NAME}_SOURCE_DIR}/pipelines"
                    ${PIPELINE_DEST_DIR}
                    COMMENT "Copying Pipeline Folder into Binary Directory"
                    COMMAND_EXPAND_LISTS
                    VERBATIM
                  )
  set_target_properties(Copy_${PLUGIN_NAME}_Pipeline_Folder PROPERTIES FOLDER Plugins/${PLUGIN_NAME})

  option(${PLUGIN_NAME}_INSTALL_PIPELINES "Enables install of ${PLUGIN_NAME} pipelines" ON)

  if(${PLUGIN_NAME}_INSTALL_PIPELINES)
    set(INSTALL_DESTINATION "pipelines/${PLUGIN_NAME}")
    install(DIRECTORY
      "${${PLUGIN_NAME}_SOURCE_DIR}/pipelines/"
      DESTINATION "${INSTALL_DESTINATION}"
      COMPONENT Applications
    )
  endif()
endif()

# -----------------------------------------------------------------------
# Download Example Data Files
# -----------------------------------------------------------------------
include(${simplnx_SOURCE_DIR}/cmake/Utility.cmake)
set(DATA_DEST_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CX_CONFIG_DIR}/Data")

file(TO_CMAKE_PATH "${DREAM3D_DATA_DIR}" DREAM3D_DATA_DIR_NORM)

# -----------------------------------------------------------------------
# These data files are pulled from the 'simplnx data archive' on github,
# decompressed and then copied to the build folder. Install rules are also
# generated for each of the data sets.
# -----------------------------------------------------------------------
if(EXISTS "${DREAM3D_DATA_DIR}" AND SIMPLNX_DOWNLOAD_TEST_FILES) 
  if(NOT EXISTS ${DREAM3D_DATA_DIR}/TestFiles/)
    file(MAKE_DIRECTORY "${DREAM3D_DATA_DIR}/TestFiles/")
  endif()
  # download_test_data(DREAM3D_DATA_DIR ${DREAM3D_DATA_DIR}
  #                  ARCHIVE_NAME T12-MAI-2010.tar.gz
  #                  SHA512 e33f224d19ad774604aa28a3263a00221a3a5909040685a3d14b6cba78e36d174b045223c28b462ab3eaea0fbc1c9f0657b1bd791a947799b9f088b13d777568
  #                  INSTALL
  #                  )
  # add_custom_target(Copy_${PLUGIN_NAME}_T12-MAI-2010 ALL
  #                  COMMAND ${CMAKE_COMMAND} -E tar xzf "${DREAM3D_DATA_DIR}/TestFiles/T12-MAI-2010.tar.gz" 
  #                  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different "${DREAM3D_DATA_DIR}/TestFiles/T12-MAI-2010" "${DATA_DEST_DIR}/T12-MAI-2010"
  #                  COMMAND ${CMAKE_COMMAND} -E rm -rf "${DREAM3D_DATA_DIR}/TestFiles/T12-MAI-2010"
  #                  WORKING_DIRECTORY "${DREAM3D_DATA_DIR}/TestFiles"
  #                  COMMENT "Copying ${PLUGIN_NAME}/T12-MAI-2010 data into Binary Directory"
  #                  DEPENDS Fetch_Remote_Data_Files  # Make sure all remote files are downloaded before trying this
  #                  COMMAND_EXPAND_LISTS
  #                  VERBATIM
  #                )
  # set_target_properties(Copy_${PLUGIN_NAME}_T12-MAI-2010 PROPERTIES FOLDER Plugins/${PLUGIN_NAME})

endif()

# -----------------------------------------------------------------------
# Create build folder copy rules and install rules for the 'data' folder
# for this plugin
# -----------------------------------------------------------------------
add_custom_target(Copy_${PLUGIN_NAME}_Data ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${${PLUGIN_NAME}_SOURCE_DIR}/data ${DATA_DEST_DIR}/${PLUGIN_NAME}
  COMMENT "Copying ${PLUGIN_NAME} data into Binary Directory"
  COMMAND_EXPAND_LISTS
  VERBATIM
)
set_target_properties(Copy_${PLUGIN_NAME}_Data PROPERTIES FOLDER Plugins/${PLUGIN_NAME})

option(${PLUGIN_NAME}_INSTALL_DATA_FILES "Enables install of ${PLUGIN_NAME} data files" ON)

set(Installed_Data_Files

)

if(${PLUGIN_NAME}_INSTALL_DATA_FILES)
  install(FILES
    ${Installed_Data_Files}
    DESTINATION Data/${PLUGIN_NAME}
    COMPONENT Applications
  )
endif()
