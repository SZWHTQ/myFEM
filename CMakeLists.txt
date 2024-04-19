cmake_minimum_required(VERSION 3.15)

project(
    myFEM
    VERSION 1.0.0
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

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

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: version ${Eigen3_VERSION}")
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

set(srcs)
add_subdirectory(src)

add_library(
    ${PROJECT_NAME}
    STATIC
    ${srcs}
)
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
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
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)
vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)
install(
    TARGETS ${PROJECT_NAME}
    DESTINATION lib
    EXPORT ${PROJECT_NAME}Targets
)

add_subdirectory(work)

enable_testing()
add_subdirectory(test)

set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

# Choose different generators, such as NSIS (Windows installer), DEB (Debian package), etc.
set(CPACK_GENERATOR "TGZ")
include(CPack)