set(SRCS
    configFile.cpp
)
add_library(configFile OBJECT ${SRCS})
target_link_libraries(configFile minivector)