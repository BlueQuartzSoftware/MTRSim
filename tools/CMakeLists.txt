add_executable(generate_pole_figures
  generate_pole_figures.cpp
)

target_include_directories(generate_pole_figures
  PRIVATE
    ${Stb_INCLUDE_DIR}
)

target_link_libraries(generate_pole_figures
  PRIVATE
    EbsdLib::EbsdLib
)
