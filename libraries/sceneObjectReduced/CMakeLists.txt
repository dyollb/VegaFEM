set(SRCS sceneObjectReduced.cpp sceneObjectReduced6DOF.cpp
         sceneObjectReducedCPU.cpp sceneObjectReducedCPU6DOF.cpp)

if(FOUND_CG)
  list(APPEND SRCS sceneObjectReducedGPU.cpp sceneObjectReducedGPU6DOF.cpp)
endif()

collect_sources(VEGA_SOURCES ${SRCS})
