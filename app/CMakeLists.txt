add_executable(mtrsim_app
  main.cpp
)

set_target_properties(mtrsim_app PROPERTIES
  OUTPUT_NAME mtrsim
)

target_link_libraries(mtrsim_app
  PRIVATE
    mtrsim::mtrsim
    CLI11::CLI11
    nlohmann_json::nlohmann_json
    spdlog::spdlog
)
