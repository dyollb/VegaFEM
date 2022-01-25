set(SRCS
    reducedForceModel.cpp
    reducedForceModelWithHessian.cpp
)
add_library(reducedForceModel OBJECT ${SRCS})
target_link_libraries(reducedForceModel PRIVATE matrixIO ${VEGAFEM_DEPENDENCIES})