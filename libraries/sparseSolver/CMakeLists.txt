set(SRCS
    # ARPACKSolver.cpp
    CGSolver.cpp
    invMKSolver.cpp
    invZTAZMSolver.cpp
    invZTAZSolver.cpp
    LagrangeMultiplierSolver.cpp
    linearSolver.cpp
    PardisoSolver.cpp
    SPOOLESSolver.cpp
    SPOOLESSolverMT.cpp
    ZTAZMultiplicator.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
