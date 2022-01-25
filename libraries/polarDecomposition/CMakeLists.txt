set(SRCS
    polarDecomposition.cpp
    polarDecompositionGradient.cpp
    polarDecompositionGradientDriver.cpp
)
add_library(polarDecomposition OBJECT ${SRCS})
target_link_libraries(polarDecomposition PRIVATE minivector)