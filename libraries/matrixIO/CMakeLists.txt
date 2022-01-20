set(SRCS
    matrixIO.cpp
    multiMatrixIO.cpp
)
add_library(matrixIO OBJECT ${SRCS})
target_link_libraries(matrixIO )