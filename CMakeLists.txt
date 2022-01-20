cmake_minimum_required(VERSION 3.1)
project(VegaFEM)

#FindOpenGL
find_package(OpenGL QUIET)
# #FindGLUT
find_package(GLUT QUIET)
# #FindGLEW
find_package(GLEW QUIET)
# #FindCg
find_package(Cg QUIET)
#FindBLAS
find_package(BLAS QUIET)
if(BLAS_FOUND)
	#FindLAPACK
	find_package(LAPACK QUIET)
endif(BLAS_FOUND)
# #FindCGAL
find_package(CGAL QUIET)
# #FindTBB
find_package(TBB QUIET)


include_directories(
	${OPENGL_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include
)

if (GLUT_FOUND)
	include_directories(${GLUT_INCLUDE_DIR})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libraries/include/openGL-headers.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libraries/include/openGL-headers.h)

if (TBB_FOUND)
	get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
	add_definitions(-DUSE_TBB)
	include_directories(${TBB_INCLUDE_DIRS})
endif()


set(VEGA_CORE_LIBRARIES
	basicAlgorithms
	clothBW
	configFile
	constrainedDOFs
	corotationalLinearFEM
	elasticForceModel
	forceModel
	getopts
	graph
	hashTable
	imageIO
	integrator
	interpolationCoordinates
	isotropicHyperelasticFEM
	listIO
	massSpringSystem
	matrixIO
	mesh
	minivector
	objMesh
	performanceCounter
	polarDecomposition
	quaternion
	rigidBodyDynamics
	sparseMatrix
	stencilForceModel
	stvk
	volumetricMesh
)

if(OPENGL_FOUND AND FOUND_CG)
	list(APPEND VEGA_CORE_LIBRARIES objMeshGPUDeformer)
endif()

if(OPENGL_FOUND AND GLUT_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		camera
		lighting
		openGLHelper
		renderVolumetricMesh
		sceneObject
	)
	if(GLEW_FOUND)
		list(APPEND VEGA_CORE_LIBRARIES glslPhong)
	endif()
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		integratorDense
		integratorSparse
		laplacianMatrix
		modalMatrix
		matrix
		reducedElasticForceModel
		reducedForceModel
		reducedStvk
		sceneObjectReduced
		sparseSolver
	)
endif()

if(CGAL_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		exactArithmetic
		virtualTets
	)
endif()


foreach(lib ${VEGA_CORE_LIBRARIES})
	include_directories(libraries/${lib})
endforeach()

add_definitions(-D_USE_MATH_DEFINES)
if (WIN32)
	add_definitions(-DNOMINMAX)
endif()

foreach(lib ${VEGA_CORE_LIBRARIES})
	add_subdirectory(libraries/${lib})
endforeach()

if(BUILD_UTILITIES)
	add_subdirectory(utilities)
endif()


# linke everything together into one library
set(OBJECT_FILES)
foreach(lib ${VEGA_CORE_LIBRARIES})
	list(APPEND OBJECT_FILES $<TARGET_OBJECTS:${lib}>)
endforeach()

add_library(vega ${OBJECT_FILES})
target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries)
target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries/include)
foreach(lib ${VEGA_CORE_LIBRARIES})
	target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib})
endforeach()