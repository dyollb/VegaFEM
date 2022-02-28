set(SRCS
    clothBWStencilForceModel.cpp
    corotationalLinearFEMStencilForceModel.cpp
    forceModelAssembler.cpp
    isotropicHyperelasticFEMStencilForceModel.cpp
    linearFEMStencilForceModel.cpp
    massSpringStencilForceModel.cpp
    StVKStencilForceModel.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
