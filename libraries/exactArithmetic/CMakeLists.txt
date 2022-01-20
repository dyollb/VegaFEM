set(SRCS
    er.cpp
    er_toDouble.cpp
    planeER.cpp
    vec3ER.cpp
)
add_library(exactArithmetic OBJECT ${SRCS})
target_link_libraries(exactArithmetic )