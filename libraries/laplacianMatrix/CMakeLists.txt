set(SRCS
    laplacianMatrix.cpp
)
add_library(laplacianMatrix OBJECT ${SRCS})
target_link_libraries(laplacianMatrix sparseMatrix volumetricMesh minivector sparseSolver)