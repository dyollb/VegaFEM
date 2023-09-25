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
    objMeshRender.cpp
    octree.cpp
    simpleSphere.cpp
    triangle-closestPoint.cpp
    triangle.cpp
    tribox3.cpp)

collect_sources(VEGA_SOURCES ${SRCS})