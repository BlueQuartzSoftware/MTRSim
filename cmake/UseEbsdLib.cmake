# ------------------------------------------------------------------------------
# Required EbsdLib and H5Support
# ------------------------------------------------------------------------------

if(MTRSIM_USE_LOCAL_EBSD_LIB)
    if(NOT TARGET EbsdLib::EbsdLib)
        if("${EbsdLibProj_SOURCE_DIR}" STREQUAL "")
            message(STATUS "EbsdLibProj_SOURCE_DIR being set to '${MTRSim_SOURCE_DIR}/../EbsdLib'")
            set(EbsdLibProj_SOURCE_DIR "${MTRSim_SOURCE_DIR}/../EbsdLib")
        else()
            message(STATUS "EbsdLibProj_SOURCE_DIR: ${EbsdLibProj_SOURCE_DIR}")
        endif()

        if(NOT EXISTS "${EbsdLibProj_SOURCE_DIR}")
            message(FATAL_ERROR "${EbsdLibProj_SOURCE_DIR} does not exist. Please set the EbsdLibProj_SOURCE_DIR variable to the EbsdLib directory.")
        endif()

        set(EbsdLib_ENABLE_HDF5 ON)
        set(EbsdLib_USE_PARALLEL_ALGORITHMS ${MTRSIM_ENABLE_MULTICORE})
        set(EbsdLib_BUILD_H5SUPPORT ON)
        set(H5Support_INCLUDE_QT_API OFF)
        add_subdirectory( ${EbsdLibProj_SOURCE_DIR} ${PROJECT_BINARY_DIR}/EbsdLib)
    endif()
else()
    find_package(H5Support REQUIRED)
    find_package(EbsdLib REQUIRED)
endif()

