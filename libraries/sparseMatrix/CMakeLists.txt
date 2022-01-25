set(SRCS
    sparseMatrix.cpp
)
add_library(sparseMatrix OBJECT ${SRCS})
target_link_libraries(sparseMatrix PRIVATE constrainedDOFs ${VEGA_DEPENDENCIES})