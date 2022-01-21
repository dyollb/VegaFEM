set(SRCS
    barycentricCoordinates.cpp
)
add_library(interpolationCoordinates OBJECT ${SRCS})
target_link_libraries(interpolationCoordinates minivector objMesh volumetricMesh mesh)