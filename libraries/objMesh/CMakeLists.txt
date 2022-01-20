set(SRCS
    createObjMesh.cpp
    objMesh.cpp
    objMeshBinaryLoader.cpp
    objMeshClose.cpp
    objMeshGraph.cpp
    objMeshOctree.cpp
    objMeshOffsetVoxels.cpp
    objMeshOrientable.cpp
    octree.cpp
)
if (GLUT_FOUND)
    list(APPEND SRCS objMeshRender.cpp)
endif()

add_library(objMesh OBJECT ${SRCS})
target_link_libraries(objMesh minivector matrixIO graph imageIO openGLHelper basicAlgorithms mesh)