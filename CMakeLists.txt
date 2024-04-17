cmake_minimum_required(VERSION 3.12)

project(myFEM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)

set(srcs)
add_subdirectory(src)

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

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toml
    ${CMAKE_SOURCE_DIR}/include
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

add_subdirectory(work)

enable_testing()
add_subdirectory(test)
