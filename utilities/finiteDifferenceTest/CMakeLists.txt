set(SRCS
    finiteDifferenceTest.cpp
    vectorHelper.cpp
)

add_executable(finiteDifferenceTest ${SRCS})
target_link_libraries(finiteDifferenceTest vega ${VEGA_DEPENDENCIES})