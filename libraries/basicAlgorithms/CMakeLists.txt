set(SRCS
    averagingBuffer.cpp
    containerHelper.cpp
    disjointSet.cpp
    fileIO.cpp
    stringHelper.cpp
    uniqueIntegerID.cpp
)
add_library(basicAlgorithms OBJECT ${SRCS})
#target_link_libraries(basicAlgorithms )