set(SRCS
    getIntegratorSolver.cpp
    integratorBase.cpp
)
add_library(integrator OBJECT ${SRCS})
target_link_libraries(integrator )