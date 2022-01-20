set(SRCS
    imageIO.cpp
)
add_library(imageIO OBJECT ${SRCS})
target_link_libraries(imageIO )