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
target_link_libraries(stencilForceModel forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM clothBW)