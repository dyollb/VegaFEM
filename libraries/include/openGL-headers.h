#include "vega-config.h"

#if defined(_WIN32) || defined(WIN32)
  #include <Windows.h>
#endif

#if defined(_WIN32) || defined(WIN32) || defined(linux) || defined (__linux__)
  #include <GL/gl.h> 
  #include <GL/glu.h>
  #if defined(VEGA_USE_GLUT)
    #include <GL/glut.h>
  #endif
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
  #if defined(VEGA_USE_GLUT)
    #include <GLUT/glut.h>
  #endif
#endif