set(SRCS
    homogeneousMooneyRivlinIsotropicMaterial.cpp
    homogeneousNeoHookeanIsotropicMaterial.cpp
    homogeneousStVKIsotropicMaterial.cpp
    isotropicHyperelasticFEM.cpp
    isotropicHyperelasticFEMMT.cpp
    isotropicMaterial.cpp
    isotropicMaterialWithCompressionResistance.cpp
    MooneyRivlinIsotropicMaterial.cpp
    neoHookeanIsotropicMaterial.cpp
    StVKIsotropicMaterial.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
