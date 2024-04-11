set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

# Check system
if(WIN32)
    set(GMSH_DIR "D:/ProgramFiles/gmsh/SDK")
else()
    set(GMSH_DIR "~/Repository/Gmsh/SDK/latest")
endif()

include_directories(
    "${GMSH_DIR}/include"
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
    "${GMSH_DIR}/lib"
)
