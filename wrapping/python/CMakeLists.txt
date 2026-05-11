pybind11_add_module(mtrsimpy mtrsim.cpp)

set_target_properties(mtrsimpy PROPERTIES
  OUTPUT_NAME mtrsim
  LIBRARY_OUTPUT_DIRECTORY $<TARGET_FILE_DIR:libMTRSim>
)

target_compile_definitions(mtrsimpy PRIVATE PYBIND11_DETAILED_ERROR_MESSAGES)

target_link_libraries(mtrsimpy PRIVATE libMTRSim::libMTRSim pybind11::headers)

target_compile_features(mtrsimpy PRIVATE cxx_std_20)

target_include_directories(mtrsimpy PRIVATE ${MTRSim_SOURCE_DIR}/src)

# CTest smoke test: verify the module imports without crashing
if(MTRSIM_BUILD_TESTS AND Python3_EXECUTABLE)
  add_test(
    NAME Python::ImportSmoke
    COMMAND ${Python3_EXECUTABLE} -c "import mtrsim; print('OK')"
  )
  set_property(TEST Python::ImportSmoke PROPERTY
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:libMTRSim>"
  )
endif()
