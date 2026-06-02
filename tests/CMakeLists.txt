add_executable(mtrsim_tests
  main.cpp
  test_config_io.cpp
  test_placeholder.cpp
  test_qsimvn.cpp
  test_gp_generator.cpp
  test_odf_file_io.cpp
  test_odf_sampler.cpp
  test_mtrsim_driver.cpp
  test_ipf_mapper.cpp
  test_odf_builder.cpp
  test_observer.cpp
)

target_include_directories(mtrsim_tests
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/LibMTRSim
)

target_link_libraries(mtrsim_tests
  PRIVATE
    libMTRSim::libMTRSim
    Catch2::Catch2
)

target_compile_definitions(mtrsim_tests PRIVATE MTRSIM_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

catch_discover_tests(mtrsim_tests)
