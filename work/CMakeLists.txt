set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Check system
if(WIN32)
    set(GMSH_DIR "D:/ProgramFiles/gmsh/SDK")
else()
    set(GMSH_DIR "~/Repository/Gmsh/SDK/latest")
endif()

add_library(
    work
    STATIC
    ${dir}/ApplyBoundary.cc
    ${dir}/GenerateMesh.cc
    ${dir}/GetStrainEnergyChange.cc
    ${dir}/SetMaterial.cc
)
target_link_libraries(
    work
    PUBLIC
    gmsh
    ${PROJECT_NAME}
)
target_link_directories(
    work
    PUBLIC
    ${GMSH_DIR}/lib
)
target_include_directories(
    work
    PUBLIC
    ${GMSH_DIR}/include
    ${dir}
)
set_target_properties(
    work PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)
install(
    TARGETS work
    DESTINATION lib
    EXPORT workTargets
)

add_library(
    pythonWorker
    SHARED
    ${dir}/worker.cc
)
target_link_libraries(
    pythonWorker
    PRIVATE
    work
)
install(
    TARGETS pythonWorker
    DESTINATION lib
    EXPORT workTargets
)

if(WIN32)
    target_compile_definitions(pythonWorker PRIVATE EXAMPLE_EXPORTS)
endif()

set_target_properties(pythonWorker PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE TRUE
    OUTPUT_NAME "worker"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/library/build
    VISIBILITY_INLINES_HIDDEN ON
    CXX_VISIBILITY_PRESET hidden
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

add_subdirectory(single)
