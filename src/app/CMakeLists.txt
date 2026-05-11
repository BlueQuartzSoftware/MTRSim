add_executable(mtrsim_app
  main.cpp
)

set_target_properties(mtrsim_app PROPERTIES
  OUTPUT_NAME libMTRSim
)

target_link_libraries(mtrsim_app
  PRIVATE
    libMTRSim::libMTRSim
    CLI11::CLI11
    HDF5::HDF5
    nlohmann_json::nlohmann_json
    spdlog::spdlog
)
