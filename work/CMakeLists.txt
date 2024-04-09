set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(
    "~/Repository/Gmsh/SDK/latest/include"
    "${CMAKE_HOME_DIRECTORY}/toml"
    "${dir}"
)

add_executable(
    main
    "${dir}/ApplyBoundary.cc"
    "${dir}/GenerateMesh.cc"
    "${dir}/GetStrainEnergyChange.cc"
    "${dir}/main.cc"
    "${dir}/SetMaterial.cc"
)

target_link_libraries(
    main
    ${PROJECT_NAME}
    gmsh
)

target_link_directories(
    main
    PUBLIC
    "~/Repository/Gmsh/SDK/latest/lib"
)
