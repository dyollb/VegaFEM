set(SRCS
    laplacianMatrix.cpp
)
add_library(laplacianMatrix OBJECT ${SRCS})
target_link_libraries(laplacianMatrix PRIVATE sparseMatrix volumetricMesh minivector sparseSolver)