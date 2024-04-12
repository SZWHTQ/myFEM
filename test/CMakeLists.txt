# Check system
if(WIN32)
    set(GMSH_DIR "D:/ProgramFiles/gmsh/SDK")
else()
    set(GMSH_DIR "~/Repository/Gmsh/SDK/latest")
endif()

include_directories(
    "${CMAKE_HOME_DIRECTORY}/work"
    "${GMSH_DIR}/include"
)

# add_executable(
#     meshTest
#     "${CMAKE_HOME_DIRECTORY}/work/GenerateMesh.cc"
#     "${CMAKE_CURRENT_SOURCE_DIR}/meshTest.cc"
# )

# target_link_libraries(
#     meshTest
#     gmsh
# )

# target_link_directories(
#     meshTest
#     PUBLIC
#     "${GMSH_DIR}/lib"
# )

# add_test(
#     NAME meshTest
#     COMMAND meshTest
# )

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