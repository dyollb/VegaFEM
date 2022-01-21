set(SRCS
    commandLineParser.cpp
    getopts.cpp
)
add_library(getopts OBJECT ${SRCS})
target_link_libraries(getopts )