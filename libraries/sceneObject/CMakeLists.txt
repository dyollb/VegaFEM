set(SRCS
    sceneGroundPlane.cpp sceneObject.cpp sceneObject6DOF.cpp
    sceneObjectDeformable.cpp sceneObjectDeformable6DOF.cpp
    sceneObjectWithRestPosition.cpp)

if(FOUND_CG)
  list(APPEND SRCS sceneObjectDeformableGPU.cpp
       sceneObjectDeformableGPU6DOF.cpp)
endif()

collect_sources(VEGA_SOURCES ${SRCS})
