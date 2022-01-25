set(SRCS
    centralDifferencesSparse.cpp
    eulerSparse.cpp
    implicitBackwardEulerSparse.cpp
    implicitNewmarkSparse.cpp
    integratorBaseSparse.cpp
)
add_library(integratorSparse OBJECT ${SRCS})
target_link_libraries(integratorSparse PRIVATE integrator performanceCounter constrainedDOFs sparseSolver forceModel)