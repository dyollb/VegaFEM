set(SRCS matrix.cpp matrixBLAS.cpp matrixLAPACK.cpp matrixPCA.cpp
         matrixProjection.cpp)
if(VEGA_USE_EXPOKIT)
  list(APPEND SRCS matrixExp.cpp)
endif()

collect_sources(VEGA_SOURCES ${SRCS})
