set(SRCS
    #ARPACKSolver.cpp
    CGSolver.cpp
    invMKSolver.cpp
    invZTAZMSolver.cpp
    invZTAZSolver.cpp
    LagrangeMultiplierSolver.cpp
    linearSolver.cpp
    PardisoSolver.cpp
    SPOOLESSolver.cpp
    SPOOLESSolverMT.cpp
    ZTAZMultiplicator.cpp
)
add_library(sparseSolver OBJECT ${SRCS})
target_link_libraries(sparseSolver PRIVATE sparseMatrix constrainedDOFs performanceCounter ${VEGA_DEPENDENCIES})