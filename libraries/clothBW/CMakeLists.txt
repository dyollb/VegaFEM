set(SRCS
    clothBW.cpp
    clothBWFromObjMesh.cpp
)
add_library(clothBW OBJECT ${SRCS})
target_link_libraries(clothBW PRIVATE minivector objMesh sparseMatrix)