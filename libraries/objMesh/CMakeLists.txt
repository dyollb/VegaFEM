set(SRCS 
    boundingBox.cpp
    objMeshBinaryLoader.cpp
    objMeshClose.cpp
    objMesh.cpp
    objMesh-disjointSet.cpp
    objMeshEncode.cpp
    objMeshGraph.cpp
    objMeshOctree.cpp
    objMeshOffsetVoxels.cpp
    objMeshOrientable.cpp
    octree.cpp
    simpleSphere.cpp
    triangle.cpp
    tribox3.cpp)

if(GLUT_FOUND)
    list(APPEND SRCS objMeshRender.cpp)
endif()
  
collect_sources(VEGA_SOURCES ${SRCS})