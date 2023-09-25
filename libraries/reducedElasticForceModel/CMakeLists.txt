set(SRCS 
    reducedLinearForceModel.cpp
    reducedLinearStVKForceModel.cpp
    reducedMassSpringSystemForceModel.cpp
    reducedMassSpringSystemForceModelWithHessian.cpp
    reducedSpringForceModel.cpp
    reducedStVKForceModel.cpp
    reducedStVKForceModelWithHessian.cpp
    reducedSubspaceStVKForceModel.cpp)

collect_sources(VEGA_SOURCES ${SRCS})