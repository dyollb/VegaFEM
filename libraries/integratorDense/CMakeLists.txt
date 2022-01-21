set(SRCS
    centralDifferencesDense.cpp
    implicitBackwardEulerDense.cpp
    implicitNewmarkDense.cpp
    implicitNewmarkDenseMulti1D.cpp
    integratorBaseDense.cpp
    integratorMulti1D.cpp
)
add_library(integratorDense OBJECT ${SRCS})
target_link_libraries(integratorDense matrix integrator performanceCounter reducedForceModel)