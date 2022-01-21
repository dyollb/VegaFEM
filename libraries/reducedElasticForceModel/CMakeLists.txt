set(SRCS
    reducedLinearForceModel.cpp
    reducedLinearStVKForceModel.cpp
    reducedMassSpringSystemForceModel.cpp
    reducedMassSpringSystemForceModelWithHessian.cpp
    reducedSpringForceModel.cpp
    reducedStVKForceModel.cpp
    reducedStVKForceModelWithHessian.cpp
    reducedSubspaceStVKForceModel.cpp
)
add_library(reducedElasticForceModel OBJECT ${SRCS})
target_link_libraries(reducedElasticForceModel reducedForceModel reducedStvk stvk massSpringSystem performanceCounter)