set(SRCS
    modalMatrix.cpp
    modalMatrixTransposed.cpp
)
add_library(modalMatrix OBJECT ${SRCS})
target_link_libraries(modalMatrix matrixIO)