set(SRCS
    iglRemeshSelfIntersection.cpp
)

add_library(libiglInterface OBJECT ${SRCS})
target_link_libraries(libiglInterface )
target_include_directories(libiglInterface PRIVATE ${IGL_INCLUDE_DIRS})