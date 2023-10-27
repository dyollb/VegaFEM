set(SRCS
    computeStiffnessMatrixNullspace.cpp
    cubicMesh.cpp
    generateGradientMatrix.cpp
    generateInterpolationMatrix.cpp
    generateMassMatrix.cpp
    generateMeshGraph.cpp
    generateSurfaceMesh.cpp
    generateTetMeshFromCubicMesh.cpp
    interpolationWeightsMultiLoad.cpp
    tetMesh.cpp
    volumetricMesh.cpp
    volumetricMeshDeformationGradient.cpp
    volumetricMeshENuMaterial.cpp
    volumetricMeshLoader.cpp
    volumetricMeshMooneyRivlinMaterial.cpp
    volumetricMeshOrthotropicMaterial.cpp
    volumetricMeshParser.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
