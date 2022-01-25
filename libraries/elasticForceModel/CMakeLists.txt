set(SRCS
    clothBWForceModel.cpp
    corotationalLinearFEMForceModel.cpp
    isotropicHyperelasticFEMForceModel.cpp
    linearFEMForceModel.cpp
    massSpringSystemForceModel.cpp
    StVKForceModel.cpp
)
add_library(elasticForceModel OBJECT ${SRCS})
target_link_libraries(elasticForceModel PRIVATE forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM clothBW)