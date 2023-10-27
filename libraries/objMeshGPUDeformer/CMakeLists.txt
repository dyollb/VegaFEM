set(SRCS
    coarseToFine-shaders.cpp
    objMeshGPUDeformer.cpp
    objMeshGPUDeformer_coarseToFine.cpp
    objMeshGPUDeformer_coarseToFine_fbo.cpp
    objMeshGPUDeformer_uUq.cpp
    objMeshGPUDeformer_uUq_fbo.cpp
    objMeshGPUDeformer_uUq_pbuffer.cpp
    objMeshGPUDeformer_uUq_setLighting.cpp
    uUq-shaders.cpp
    vbo.cpp)

collect_sources(VEGA_SOURCES ${SRCS})
