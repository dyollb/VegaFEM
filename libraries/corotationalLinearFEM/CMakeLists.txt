set(SRCS
    corotationalLinearFEM.cpp
)
add_library(corotationalLinearFEM OBJECT ${SRCS})
target_link_libraries(corotationalLinearFEM polarDecomposition volumetricMesh sparseMatrix)