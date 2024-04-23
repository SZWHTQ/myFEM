set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Check system
if(WIN32)
    set(GMSH_DIR "$ENV{VCPKG_ROOT}/packages/gmsh_arm64-windows")
    # x64 or arm64
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "AMD64")
        set(GMSH_DIR "$ENV{VCPKG_ROOT}/packages/gmsh_x64-windows")
    elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "ARM64")
        set(GMSH_DIR "$ENV{VCPKG_ROOT}/packages/gmsh_arm64-windows")
    else()
        message(FATAL_ERROR "Unknown processor")
    endif()
elseif(LINUX OR APPLE)
    set(GMSH_DIR "~/Repository/Gmsh/SDK/latest")
else()
    message(FATAL_ERROR "Unknown system")
endif()

add_library(
    work
    SHARED
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
    PyWorker
    SHARED
    ${dir}/PyWorker.cc
)
target_link_libraries(
    PyWorker
    PRIVATE
    work
)
install(
    TARGETS PyWorker
    DESTINATION lib
    EXPORT workTargets
)

if(WIN32)
    target_compile_definitions(PyWorker PRIVATE EXAMPLE_EXPORTS)
endif()

set_target_properties(PyWorker PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE TRUE
    OUTPUT_NAME "PythonWorker"
    # LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/library/build
    VISIBILITY_INLINES_HIDDEN ON
    CXX_VISIBILITY_PRESET hidden
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

add_subdirectory(singleton)
