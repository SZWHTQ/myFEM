set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(
    "${CMAKE_HOME_DIRECTORY}/toml"
)

add_executable(
    main
    "${dir}/GetStrainEnergyChange.cc"
    "${dir}/main.cc"
)

target_link_libraries(
    main
    work
)
