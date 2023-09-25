set(SRCS
    massSpringSystem.cpp
    massSpringSystemFromCubicMesh.cpp
    massSpringSystemFromCubicMeshConfigFile.cpp
    massSpringSystemFromObjMesh.cpp
    massSpringSystemFromObjMeshConfigFile.cpp
    massSpringSystemFromTetMesh.cpp
    massSpringSystemFromTetMeshConfigFile.cpp)

if(GLUT_FOUND)
    list(APPEND SRCS renderSprings.cpp)
endif()
  

collect_sources(VEGA_SOURCES ${SRCS})