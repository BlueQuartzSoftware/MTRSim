set(MTRSIM_HEADERS
  AssignmentRule.hpp
  CrystalSymmetry.hpp
  GPGenerator.hpp
  IPFMapper.hpp
  MTRDataLoader.hpp
  ODFCalculator.hpp
  ODFSampler.hpp
  PGRFSimulation.hpp
  PoleFigure.hpp
  QSimVN.hpp
)

set(MTRSIM_SOURCES
  AssignmentRule.cpp
  CrystalSymmetry.cpp
  GPGenerator.cpp
  IPFMapper.cpp
  MTRDataLoader.cpp
  ODFCalculator.cpp
  ODFSampler.cpp
  PGRFSimulation.cpp
  PoleFigure.cpp
  QSimVN.cpp
)


add_library(mtrsim SHARED ${MTRSIM_HEADERS} ${MTRSIM_SOURCES})

add_library(mtrsim::mtrsim ALIAS mtrsim)

include(GenerateExportHeader)
generate_export_header(mtrsim)

target_include_directories(mtrsim
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/mtrsim>
  PRIVATE
    ${Stb_INCLUDE_DIR}
)

target_link_libraries(mtrsim
  PUBLIC
    Eigen3::Eigen
    HDF5::HDF5
    EbsdLib::EbsdLib
    spdlog::spdlog
  PRIVATE
    nlohmann_json::nlohmann_json
)

target_compile_features(mtrsim PUBLIC cxx_std_20)

if(MTRSIM_ENABLE_MULTICORE)
  target_link_libraries(mtrsim PUBLIC TBB::tbb)
  target_compile_definitions(mtrsim PUBLIC MTRSIM_HAS_TBB=1)
endif()
