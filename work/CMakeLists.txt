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

add_library(
    work
    SHARED
    "${dir}/ApplyBoundary.cc"
    "${dir}/GenerateMesh.cc"
    "${dir}/GetStrainEnergyChange.cc"
    "${dir}/SetMaterial.cc"
)
target_link_libraries(
    work
    gmsh
    ${PROJECT_NAME}
)

target_link_directories(
    work
    PUBLIC
    "${GMSH_DIR}/lib"
)

add_library(
    Worker
    SHARED
    "${dir}/worker.cc"
)
target_link_libraries(
    Worker
    work
)

if(WIN32)
    target_compile_definitions(Worker PRIVATE "EXAMPLE_EXPORTS")
endif()

set_target_properties(Worker PROPERTIES
    PREFIX ""
    OUTPUT_NAME "worker"
    POSITION_INDEPENDENT_CODE TRUE
    VISIBILITY_INLINES_HIDDEN ON
    CXX_VISIBILITY_PRESET hidden
)
target_compile_options(Worker PRIVATE "-fdeclspec")

add_subdirectory("single")
