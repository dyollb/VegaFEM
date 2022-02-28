set(SRCS
    massSpringSystem.cpp
    massSpringSystemFromCubicMesh.cpp
    massSpringSystemFromCubicMeshConfigFile.cpp
    massSpringSystemFromObjMesh.cpp
    massSpringSystemFromObjMeshConfigFile.cpp
    massSpringSystemFromTetMesh.cpp
    massSpringSystemFromTetMeshConfigFile.cpp
    renderSprings.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
