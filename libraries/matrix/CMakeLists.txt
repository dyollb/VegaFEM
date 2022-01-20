set(SRCS
    matrix.cpp
    matrixBLAS.cpp
    matrixBLASOptimized.cpp
    matrixBLASVanilla.cpp
    matrixExp.cpp
    matrixLAPACK.cpp
    matrixPCA.cpp
    matrixProjection.cpp
)
add_library(matrix OBJECT ${SRCS})
target_link_libraries(matrix matrixIO)