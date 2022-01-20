set(SRCS
    eig3.cpp
    mat3d.cpp
    vec2d.cpp
    vec3.cpp
    vec3d.cpp
    vec3i.cpp
    vec4i.cpp
)
add_library(minivector OBJECT ${SRCS})
target_link_libraries(minivector )