set(MTRSIM_HEADERS
  AssignmentRule.hpp
  ConfigIO.hpp
  GPGenerator.hpp
  IPFMapper.hpp
  ISimulationObserver.hpp
  MTRDataLoader.hpp
  MTRSimDriver.hpp
  ODFBuilder.hpp
  ODFCalculator.hpp
  ODFFileIO.hpp
  ODFSampler.hpp
  PGRFSimulation.hpp
  PoleFigure.hpp
  QSimVN.hpp
  SimulationObservers.hpp
)

set(MTRSIM_SOURCES
  AssignmentRule.cpp
  ConfigIO.cpp
  GPGenerator.cpp
  IPFMapper.cpp
  MTRDataLoader.cpp
  MTRSimDriver.cpp
  ODFBuilder.cpp
  ODFCalculator.cpp
  ODFFileIO.cpp
  ODFSampler.cpp
  PGRFSimulation.cpp
  PoleFigure.cpp
  QSimVN.cpp
)


add_library(libMTRSim SHARED ${MTRSIM_HEADERS} ${MTRSIM_SOURCES})

add_library(libMTRSim::libMTRSim ALIAS libMTRSim)

include(GenerateExportHeader)
generate_export_header(libMTRSim)

target_include_directories(libMTRSim
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>

    $<INSTALL_INTERFACE:include/libMTRSim>
  PRIVATE
    ${Stb_INCLUDE_DIR}
)

target_link_libraries(libMTRSim
  PUBLIC
    Eigen3::Eigen
    HDF5::HDF5
    EbsdLib::EbsdLib
    spdlog::spdlog
  PRIVATE
    nlohmann_json::nlohmann_json
)

target_compile_features(libMTRSim PUBLIC cxx_std_20)

if(MTRSIM_ENABLE_MULTICORE)
  target_link_libraries(libMTRSim PUBLIC TBB::tbb)
  target_compile_definitions(libMTRSim PUBLIC MTRSIM_HAS_TBB=1)
endif()
