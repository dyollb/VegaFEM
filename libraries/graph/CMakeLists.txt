set(SRCS
    graph.cpp
)
add_library(graph OBJECT ${SRCS})
target_link_libraries(graph matrixIO sparseMatrix)