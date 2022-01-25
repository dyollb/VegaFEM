set(SRCS
    clothBWStencilForceModel.cpp
    corotationalLinearFEMStencilForceModel.cpp
    forceModelAssembler.cpp
    isotropicHyperelasticFEMStencilForceModel.cpp
    linearFEMStencilForceModel.cpp
    massSpringStencilForceModel.cpp
    StVKStencilForceModel.cpp
)
add_library(stencilForceModel OBJECT ${SRCS})
target_link_libraries(stencilForceModel PRIVATE forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM clothBW ${VEGA_DEPENDENCIES})