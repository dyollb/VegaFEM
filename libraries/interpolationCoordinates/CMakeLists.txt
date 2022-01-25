set(SRCS
    barycentricCoordinates.cpp
)
add_library(interpolationCoordinates OBJECT ${SRCS})
target_link_libraries(interpolationCoordinates PRIVATE minivector objMesh volumetricMesh mesh)