set(SRCS
    constrainedDOFs.cpp
)
add_library(constrainedDOFs OBJECT ${SRCS})
target_link_libraries(constrainedDOFs )