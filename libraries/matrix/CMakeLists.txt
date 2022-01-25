set(SRCS
    matrix.cpp
    matrixBLAS.cpp
    matrixExp.cpp
    matrixLAPACK.cpp
    matrixPCA.cpp
    matrixProjection.cpp
)
add_library(matrix OBJECT ${SRCS})
target_link_libraries(matrix PRIVATE matrixIO ${VEGAFEM_DEPENDENCIES})