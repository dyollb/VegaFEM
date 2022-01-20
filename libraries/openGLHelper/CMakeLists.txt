set(SRCS
    handleRender.cpp
    inputDevice.cpp
    openGLHelper.cpp
    printBitmap.cpp
    saveScreenShot.cpp
)
add_library(openGLHelper OBJECT ${SRCS})
target_link_libraries(openGLHelper imageIO minivector performanceCounter)