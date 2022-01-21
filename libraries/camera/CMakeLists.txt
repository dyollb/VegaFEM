set(SRCS
    camera.cpp
    cameraChangeLoad.cpp
)
add_library(camera OBJECT ${SRCS})
target_link_libraries(camera )