set(SRCS
    hashTable.cpp
)
add_library(hashTable OBJECT ${SRCS})
target_link_libraries(hashTable )