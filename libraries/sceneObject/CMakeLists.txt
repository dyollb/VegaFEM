set(SRCS
    sceneGroundPlane.cpp
    sceneObject.cpp
    sceneObject6DOF.cpp
    sceneObjectDeformable.cpp
    sceneObjectDeformable6DOF.cpp
    sceneObjectWithRestPosition.cpp
)

if (FOUND_CG)
    list(APPEND SRCS
        sceneObjectDeformableGPU.cpp
        sceneObjectDeformableGPU6DOF.cpp
    )
endif()

add_library(sceneObject OBJECT ${SRCS})
target_link_libraries(sceneObject objMesh lighting minivector imageIO)