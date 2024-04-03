include_directories(
    "${CMAKE_HOME_DIRECTORY}/work"
    "~/Repository/Gmsh/SDK/latest/include"
)

add_executable(
    meshTest
    "${CMAKE_HOME_DIRECTORY}/work/GenerateMesh.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/meshTest.cc"
)

target_link_libraries(
    meshTest
    gmsh
)

target_link_directories(
    meshTest
    PUBLIC
    "~/Repository/Gmsh/SDK/latest/lib"
)

add_test(
    NAME meshTest
    COMMAND meshTest
)

add_executable(
    elementKTest
    "${CMAKE_HOME_DIRECTORY}/test/elementStiffnessMatrix.cc"
)

target_link_libraries(
    elementKTest
    ${PROJECT_NAME}
)

add_test(
    NAME elementKTest
    COMMAND elementKTest
)