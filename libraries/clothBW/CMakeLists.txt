set(SRCS
    clothBW.cpp
    clothBWFromObjMesh.cpp
)
add_library(clothBW OBJECT ${SRCS})
target_link_libraries(clothBW minivector objMesh sparseMatrix)