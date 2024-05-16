set(dir ${CMAKE_CURRENT_SOURCE_DIR})

list(
    APPEND srcs
    ${dir}/Element.cc
    ${dir}/Mesh.cc
    ${dir}/Serendipity.cc
    ${dir}/ThreadPool.cc
    ${dir}/Timer.cc
    ${dir}/vtkManager.cc
)

set(srcs ${srcs} PARENT_SCOPE)