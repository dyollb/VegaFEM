set(SRCS
    quaternion.cpp
)
add_library(quaternion OBJECT ${SRCS})
target_link_libraries(quaternion )