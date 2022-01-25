set(SRCS
    profiler.cpp
)
add_library(performanceCounter OBJECT ${SRCS})
target_link_libraries(performanceCounter PRIVATE matrixIO basicAlgorithms)