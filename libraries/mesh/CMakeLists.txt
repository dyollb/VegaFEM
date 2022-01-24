set(SRCS
    boundingBox.cpp
    createTriMesh.cpp
    exactOctree.cpp
    geometryQuery.cpp
    intersection_tunicate.cpp
    labelOuterTets.cpp
    meshIntersection.cpp
    predicates.cpp
    rectKey.cpp
    simpleSphere.cpp
    tetKey.cpp
    tetMeshGeo.cpp
    tetMeshManifold.cpp
    triangle.cpp
    tribox3.cpp
    triKey.cpp
    triMeshGeo.cpp
    triMeshManifold.cpp
    triMeshNeighbor.cpp
    triMeshPseudoNormal.cpp
    verticesInfo.cpp
    windingNumberTree.cpp
)
add_library(mesh OBJECT ${SRCS})
target_link_libraries(mesh minivector basicAlgorithms ${VEGAFEM_DEPENDENCIES})