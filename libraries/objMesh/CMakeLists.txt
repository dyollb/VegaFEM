set(SRCS
    createObjMesh.cpp
    objMesh.cpp
    objMeshBinaryLoader.cpp
    objMeshClose.cpp
    objMeshGraph.cpp
    objMeshOctree.cpp
    objMeshOffsetVoxels.cpp
    objMeshOrientable.cpp
    octree.cpp)
if(GLUT_FOUND)
  list(APPEND SRCS objMeshRender.cpp)
endif()

collect_sources(VEGA_SOURCES ${SRCS})
