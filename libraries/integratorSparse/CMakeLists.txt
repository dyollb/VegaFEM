set(SRCS
    centralDifferencesSparse.cpp eulerSparse.cpp
    implicitBackwardEulerSparse.cpp implicitNewmarkSparse.cpp
    integratorBaseSparse.cpp)

collect_sources(VEGA_SOURCES ${SRCS})