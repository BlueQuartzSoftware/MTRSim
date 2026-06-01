# ------------------------------------------------------------------------------
# If there are additional library packages that need to be found, this is where
# to do that using the usual find_package(....) call
# find_package(....)

# ------------------------------------------------------------------------------
# Define the list of unit test source files
set(${PLUGIN_NAME}UnitTest_SRCS
    ComputeODFTest.cpp
    MTRSimTest.cpp
    ReadMTRSimODFTest.cpp
    WriteMTRSimODFTest.cpp
)

set(DISABLED_TESTS

)

create_simplnx_plugin_unit_test(PLUGIN_NAME ${PLUGIN_NAME}
  FILTER_LIST ${${PLUGIN_NAME}UnitTest_SRCS})

# ------------------------------------------------------------------------------
# If there are additional libraries that this plugin needs to link against you
# can use the target_link_libraries() cmake call
# target_link_libraries(${PLUGIN_NAME}UnitTest PUBLIC [name of library])

# Realistic-data validation test (ComputeODFTest.cpp's [validation]-tagged
# test) loads data/real_world_microtexture_data.dream3d through SimplnxCore's
# ReadDREAM3DFilter, so we need its headers + link.
target_link_libraries(${PLUGIN_NAME}UnitTest PRIVATE SimplnxCore)
target_include_directories(${PLUGIN_NAME}UnitTest PRIVATE
  ${simplnx_SOURCE_DIR}/src/Plugins/SimplnxCore/src
)

# ------------------------------------------------------------------------------
# If there are additional source files that need to be compiled for this plugin
# can use the target_sources(...) cmake call
target_sources(${PLUGIN_NAME}UnitTest

  # PUBLIC
  # additional header files
  PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/${PLUGIN_NAME}TestUtils.hpp
)

# ------------------------------------------------------------------------------
# If there are additional include directories that are needed for this plugin
# you can use the target_include_directories(.....) cmake call
# target_include_directories(${PLUGIN_NAME}UnitTest
# PUBLIC
# ${CMAKE_CURRENT_LIST_DIR}/../../src
# )

# -----------------------------------------------------------------------------
# Here we are going to setup to download and decompress the test files. In order
# to add your own test files you will need to tar.gz the test file, compute the
# SHA 512 Hash of the file and then upload the file to
# https://github.com/bluequartzsoftware/simplnx/releases/tag/Data_Archive.
#
# Go to the web site above, "edit" the release, add the filename and SHA 512 to
# the table and then upload your compressed file.
# Save the release so that the repo is updated
# -----------------------------------------------------------------------------
option(SIMPLNX_DOWNLOAD_TEST_FILES "Download the test files" ON)

if(EXISTS "${DREAM3D_DATA_DIR}" AND SIMPLNX_DOWNLOAD_TEST_FILES)
  if(NOT EXISTS ${DREAM3D_DATA_DIR}/TestFiles/)
    file(MAKE_DIRECTORY "${DREAM3D_DATA_DIR}/TestFiles/")
  endif()

#  download_test_data(DREAM3D_DATA_DIR ${DREAM3D_DATA_DIR} ARCHIVE_NAME 12_IN625_GBCD.tar.gz SHA512 f696a8af181505947e6fecfdb1a11fda6c762bba5e85fea8d484b1af00bf18643e1d930d48f092ee238d1c19c9ce7c4fb5a8092d17774bda867961a1400e9cea)

endif()

# -----------------------------------------------------------------------------
# Include any example pipelines into the unit testing
# -----------------------------------------------------------------------------
include(${simplnx_SOURCE_DIR}/cmake/Utility.cmake)

set(TEST_PIPELINE_LIST_FILE ${${PLUGIN_NAME}_BINARY_DIR}/PrebuiltPipelineTests/nxrunner_tests.txt)
FILE(WRITE ${TEST_PIPELINE_LIST_FILE} "Prebuilt Pipeline Tests\n")
get_target_property(PIPELINE_RUNNER_NAME nxrunner NAME)
get_target_property(PIPELINE_RUNNER_DEBUG nxrunner DEBUG_POSTFIX)

