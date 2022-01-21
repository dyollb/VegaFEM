set(SRCS
    cameraLighting.cpp
    lighting.cpp
)
add_library(lighting OBJECT ${SRCS})
target_link_libraries(lighting configFile camera)