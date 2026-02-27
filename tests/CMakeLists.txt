add_executable(mtrsim_tests
  test_placeholder.cpp
)

target_link_libraries(mtrsim_tests
  PRIVATE
    mtrsim::mtrsim
    Catch2::Catch2WithMain
)

catch_discover_tests(mtrsim_tests)
