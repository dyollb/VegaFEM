set(SRCS
    massSpringSystem.cpp
    massSpringSystemFromCubicMesh.cpp
    massSpringSystemFromCubicMeshConfigFile.cpp
    massSpringSystemFromObjMesh.cpp
    massSpringSystemFromObjMeshConfigFile.cpp
    massSpringSystemFromTetMesh.cpp
    massSpringSystemFromTetMeshConfigFile.cpp
    renderSprings.cpp
)
add_library(massSpringSystem OBJECT ${SRCS})
target_link_libraries(massSpringSystem PRIVATE objMesh volumetricMesh configFile sparseMatrix)