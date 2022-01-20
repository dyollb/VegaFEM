set(SRCS
    sceneObjectReduced.cpp
    sceneObjectReduced6DOF.cpp
    sceneObjectReducedCPU.cpp
    sceneObjectReducedCPU6DOF.cpp
    sceneObjectReducedGPU.cpp
    sceneObjectReducedGPU6DOF.cpp
)
add_library(sceneObjectReduced OBJECT ${SRCS})
target_link_libraries(sceneObjectReduced objMesh lighting minivector modalMatrix sceneObject)