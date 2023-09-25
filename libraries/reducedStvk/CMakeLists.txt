set(SRCS 
    StVKReducedHessianTensor.cpp
    StVKReducedInternalForces.cpp
    StVKReducedInternalForcesMT.cpp
    StVKReducedInternalForcesMultiLoad.cpp
    StVKReducedStiffnessMatrix.cpp)

collect_sources(VEGA_SOURCES ${SRCS})