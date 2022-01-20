set(SRCS
    eig3.cpp
    mat3d.cpp
)
add_library(minivector OBJECT ${SRCS})
target_link_libraries(minivector )