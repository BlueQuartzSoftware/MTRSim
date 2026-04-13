add_executable(mtrsim_tests
  main.cpp
  test_placeholder.cpp
  test_qsimvn.cpp
  test_gp_generator.cpp
  test_odf_sampler.cpp
  test_ipf_mapper.cpp
)

target_include_directories(mtrsim_tests
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/libmtrsim
)

target_link_libraries(mtrsim_tests
  PRIVATE
    mtrsim::mtrsim
    Catch2::Catch2
)

catch_discover_tests(mtrsim_tests)
