set(SRCS
    StVKReducedHessianTensor.cpp
    StVKReducedInternalForces.cpp
    StVKReducedInternalForcesMultiLoad.cpp
    StVKReducedStiffnessMatrix.cpp
)
add_library(reducedStvk OBJECT ${SRCS})
target_link_libraries(reducedStvk PRIVATE stvk matrix ${VEGAFEM_DEPENDENCIES})