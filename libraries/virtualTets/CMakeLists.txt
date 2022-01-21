set(SRCS
    geometryQueryER.cpp
    tetTriCutFeature.cpp
    tetTriMeshCutting.cpp
    tetTriPiece.cpp
    virtualTets-via-csg.cpp
    virtualTets.cpp
)
add_library(virtualTets OBJECT ${SRCS})
target_link_libraries(virtualTets interpolationCoordinates exactArithmetic mesh basicAlgorithms libiglInterface)