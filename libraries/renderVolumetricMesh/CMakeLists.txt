set(SRCS
    renderVolumetricMesh.cpp
)
add_library(renderVolumetricMesh OBJECT ${SRCS})
target_link_libraries(renderVolumetricMesh volumetricMesh openGLHelper)