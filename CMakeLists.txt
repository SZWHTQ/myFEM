cmake_minimum_required(VERSION 3.15)

project(
    myFEM
    VERSION 1.0.0
    DESCRIPTION "A simple FEM library"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
else()
    message(WARNING "Unknown compiler, no warning flags set")
    message(WARNING "Compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Intel
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(MKLROOT $ENV{MKLROOT})
    add_definitions(-DEIGEN_USE_MKL_ALL)
    include_directories(${MKLROOT}/include)
    set(CMAKE_LINK_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl")
endif()

find_package(VTK COMPONENTS
    CommonCore
    IOCore
    IOLegacy
    IOXML
    REQUIRED
)
message(STATUS "VTK found: version ${VTK_VERSION}")

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen3 found: version ${Eigen3_VERSION}")

find_package(OpenMP)

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

if(OpenMP_CXX_FOUND)
    target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE
        OpenMP::OpenMP_CXX
    )
endif()

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