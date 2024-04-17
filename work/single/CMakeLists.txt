set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(
    main
    ${dir}/GetStrainEnergyChange.cc
    ${dir}/main.cc
)

target_link_libraries(
    main
    PRIVATE
    work
)

set_target_properties(
    main
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)
