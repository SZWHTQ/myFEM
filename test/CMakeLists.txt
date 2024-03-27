include_directories(
    "../work"
    "~/Repository/Gmsh/SDK/latest/include"
)

add_executable(
    meshTest
    "../work/GenerateMesh.cc"
    "meshTest.cc"
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