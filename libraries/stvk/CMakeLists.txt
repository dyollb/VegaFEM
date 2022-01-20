set(SRCS
    cubicMeshIntegrals.cpp
    StVKCubeABCD.cpp
    StVKElementABCD.cpp
    StVKElementABCDLoader.cpp
    StVKFEM.cpp
    StVKHessianTensor.cpp
    StVKInternalForces.cpp
    StVKStiffnessMatrix.cpp
    StVKTetABCD.cpp
    StVKTetHighMemoryABCD.cpp
)
add_library(stvk OBJECT ${SRCS})
target_link_libraries(stvk minivector volumetricMesh sparseMatrix)