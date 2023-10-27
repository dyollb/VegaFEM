set(SRCS
    centralDifferencesDense.cpp implicitBackwardEulerDense.cpp
    implicitNewmarkDense.cpp implicitNewmarkDenseMulti1D.cpp
    integratorBaseDense.cpp integratorMulti1D.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
