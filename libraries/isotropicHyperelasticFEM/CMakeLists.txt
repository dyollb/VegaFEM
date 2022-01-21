set(SRCS
    homogeneousMooneyRivlinIsotropicMaterial.cpp
    homogeneousNeoHookeanIsotropicMaterial.cpp
    homogeneousStVKIsotropicMaterial.cpp
    isotropicHyperelasticFEM.cpp
    isotropicMaterial.cpp
    isotropicMaterialWithCompressionResistance.cpp
    MooneyRivlinIsotropicMaterial.cpp
    neoHookeanIsotropicMaterial.cpp
    StVKIsotropicMaterial.cpp
)
add_library(isotropicHyperelasticFEM OBJECT ${SRCS})
target_link_libraries(isotropicHyperelasticFEM minivector volumetricMesh sparseMatrix)