include_directories(
    "~/Repository/Gmsh/SDK/latest/include"
    "./"
)

add_executable(
    main
    "./ApplyBoundary.cc"
    "./GenerateMesh.cc"
    "./main.cc"
    "./SetMaterial.cc"
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
