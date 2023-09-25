set(SRCS
    clothBWForceModel.cpp corotationalLinearFEMForceModel.cpp
    isotropicHyperelasticFEMForceModel.cpp linearFEMForceModel.cpp
    massSpringSystemForceModel.cpp StVKForceModel.cpp)

collect_sources(VEGA_SOURCES ${SRCS})