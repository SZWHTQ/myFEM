# Check system
if(WIN32)
    set(GMSH_DIR "$ENV{VCPKG_DIR}/packages/gmsh_arm64-windows")
else()
    set(GMSH_DIR "~/Repository/Gmsh/SDK/latest")
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