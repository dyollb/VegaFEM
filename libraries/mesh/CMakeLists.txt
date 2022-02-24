set(SRCS
    boundingBox.cpp
    createTriMesh.cpp
    exactOctree.cpp
    geometryQuery.cpp
    initPredicates.cpp
    intersection_tunicate.cpp
    labelOuterTets.cpp
    meshIntersection.cpp
    predicates.cpp
    predicates_Shewchuk.c
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
target_link_libraries(mesh PRIVATE minivector basicAlgorithms ${VEGA_DEPENDENCIES})