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

## Add test
# Add test for the mesh generation
add_executable(
    meshGenerationTest
    ${CMAKE_SOURCE_DIR}/work/GenerateMesh.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/meshTest.cc
)
target_link_libraries(
    meshGenerationTest
    PRIVATE
    gmsh
    ${PROJECT_NAME}
)
target_link_directories(
    meshGenerationTest
    PRIVATE
    ${GMSH_DIR}/lib
)
target_include_directories(
    meshGenerationTest
    PRIVATE
    ${GMSH_DIR}/include
    ${CMAKE_SOURCE_DIR}/work
)
add_test(
    NAME meshGeneration
    COMMAND meshGenerationTest
)

# Add test for the element stiffness matrix
add_executable(
    elementStiffnessMatrixTest
    ${CMAKE_SOURCE_DIR}/test/elementStiffnessMatrix.cc
)
target_link_libraries(
    elementStiffnessMatrixTest
    ${PROJECT_NAME}
)
add_test(
    NAME elementStiffnessMatrix
    COMMAND elementStiffnessMatrixTest
)