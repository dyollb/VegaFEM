set(SRCS
    homogeneousMooneyRivlinIsotropicMaterial.cpp
    homogeneousNeoHookeanIsotropicMaterial.cpp
    homogeneousStVKIsotropicMaterial.cpp
    isotropicHyperelasticFEM.cpp
    isotropicMaterial.cpp
    isotropicMaterialWithCompressionResistance.cpp
    MooneyRivlinIsotropicMaterial.cpp
    neoHookeanIsotropicMaterial.cpp
    StVKIsotropicMaterial.cpp)

collect_sources(VEGA_SOURCES ${SRCS})