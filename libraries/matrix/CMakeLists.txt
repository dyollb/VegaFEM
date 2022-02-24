set(SRCS
    matrix.cpp
    matrixBLAS.cpp
    matrixLAPACK.cpp
    matrixPCA.cpp
    matrixProjection.cpp
)
if (VEGA_USE_EXPOKIT)
    list(APPEND SRCS matrixExp.cpp)
endif()

add_library(matrix OBJECT ${SRCS})
target_link_libraries(matrix PRIVATE matrixIO ${VEGA_DEPENDENCIES})