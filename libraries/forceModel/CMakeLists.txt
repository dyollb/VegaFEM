set(SRCS
    finiteDifferenceTester.cpp
    forceModel.cpp
)
add_library(forceModel OBJECT ${SRCS})
target_link_libraries(forceModel sparseMatrix)