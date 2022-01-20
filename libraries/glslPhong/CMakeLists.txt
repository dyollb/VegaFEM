set(SRCS
    glslPhong.cpp
)
add_library(glslPhong OBJECT ${SRCS})
target_link_libraries(glslPhong )