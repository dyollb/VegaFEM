set(SRCS 
    cubicMeshIntegrals.cpp
    StVKCubeABCD.cpp
    StVKElementABCD.cpp
    StVKElementABCDLoader.cpp
    StVKHessianTensor.cpp
    StVKInternalForces.cpp
    StVKInternalForcesMT.cpp
    StVKStiffnessMatrix.cpp
    StVKStiffnessMatrixMT.cpp
    StVKTetABCD.cpp
    StVKTetHighMemoryABCD.cpp)

collect_sources(VEGA_SOURCES ${SRCS})