set(SRCS
    finiteDifferenceTester.cpp
    forceModel.cpp
)
add_library(forceModel OBJECT ${SRCS})
target_link_libraries(forceModel PRIVATE sparseMatrix ${VEGAFEM_DEPENDENCIES})