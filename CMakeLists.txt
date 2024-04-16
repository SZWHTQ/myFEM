cmake_minimum_required(VERSION 3.12)

project(myFEM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "~/Repository/VTK/latest/lib/cmake/vtk-9.3")
# set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "/opt/homebrew/share/eigen3/cmake")

set(srcs)
add_subdirectory("src")
include_directories(
    "include"
)

find_package(VTK COMPONENTS
    CommonCore
    IOCore
    IOLegacy
    IOXML
    REQUIRED
)
if(VTK_FOUND)
    message(STATUS "VTK found: version ${VTK_VERSION}")
else()
    message(FATAL_ERROR "VTK not found")
endif()

find_package(Eigen3)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: version ${Eigen3_VERSION}")
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

add_library(
    ${PROJECT_NAME}
    ${srcs}
)

target_link_libraries(
    ${PROJECT_NAME}
    ${VTK_LIBRARIES}
)

add_subdirectory("work")

add_subdirectory("test")
