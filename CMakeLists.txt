cmake_minimum_required(VERSION 3.18)
project(VegaFEM CXX)
set(CMAKE_CXX_STANDARD 14)


# configuration options
option(VEGA_USE_TBB "Build VegaFEM using TBB" ON)
option(VEGA_USE_GLEW "Build OpenGL shaders using GLEW" OFF)
option(VEGA_USE_GLUT "Build components that depend on GLUT" OFF)
option(VEGA_USE_CG_TOOLKIT "Build GPU components using Cg Toolkit" OFF)
option(VEGA_USE_CGAL "Build meshing components using CGAL" OFF)
option(VEGA_BUILD_UTILITIES "Build utilities" OFF)


# thirdparties
set(VEGA_DEPENDENCIES)

find_package(OpenGL REQUIRED)
list(APPEND VEGA_DEPENDENCIES OpenGL::GL)

if(VEGA_USE_TBB)
	find_package(TBB REQUIRED)
	list(APPEND VEGA_DEPENDENCIES ${TBB_IMPORTED_TARGETS})
	add_definitions(-DUSE_TBB)
endif()

if(VEGA_USE_GLUT)
	find_package(GLUT REQUIRED)
	list(APPEND VEGA_DEPENDENCIES GLUT::GLUT)
endif()

if(VEGA_USE_GLEW)
	find_package(GLEW REQUIRED)
	list(APPEND VEGA_DEPENDENCIES GLEW::GLEW) # glew_s for static
endif()

if(VEGA_USE_CG_TOOLKIT)
	find_package(Cg QUIET) # todo: this does not work yet
endif()

if(VEGA_USE_CGAL)
	find_package(CGAL REQUIRED)
	find_path(IGL_INCLUDE_DIRS igl/harmonic.h PATHS $ENV{IGL_ROOT}/include REQUIRED)
endif()

set(ENABLE_BLAS95 ON)
set(ENABLE_LAPACK95 ON)
find_package(MKL QUIET)
if(MKL_FOUND)
	set(VEGA_USE_INTEL_MKL ON)
	list(APPEND VEGA_DEPENDENCIES MKL::MKL)
else()
	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)
	list(APPEND VEGA_DEPENDENCIES BLAS::BLAS LAPACK::LAPACK)
endif()


configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include/vega-config.h.in 
	${CMAKE_BINARY_DIR}/vega-config.h
)


set(VEGA_LIBRARIES
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
	integratorDense
	integratorSparse
	interpolationCoordinates
	isotropicHyperelasticFEM
	laplacianMatrix
	listIO
	matrix
	massSpringSystem
	matrixIO
	mesh
	minivector
	modalMatrix
	objMesh
	performanceCounter
	polarDecomposition
	quaternion
	reducedElasticForceModel
	reducedForceModel
	reducedStvk
	rigidBodyDynamics
	sparseMatrix
	sparseSolver
	stencilForceModel
	stvk
	volumetricMesh
)

if(VEGA_USE_GLUT)
	list(APPEND VEGA_LIBRARIES
		camera
		lighting
		openGLHelper
		renderVolumetricMesh
		sceneObject
		sceneObjectReduced
	)
endif()

if(VEGA_USE_GLUT AND VEGA_USE_GLEW)
	list(APPEND VEGA_LIBRARIES 
		glslPhong
	)
endif()

if(VEGA_USE_CGAL)
	list(APPEND VEGA_LIBRARIES
		exactArithmetic
		libiglInterface
		virtualTets
	)
endif()

if(VEGA_USE_CG_TOOLKIT)
	list(APPEND VEGA_LIBRARIES objMeshGPUDeformer)
endif()


# include directories
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include
	${CMAKE_BINARY_DIR}
)

foreach(lib ${VEGA_LIBRARIES})
	include_directories(libraries/${lib})
endforeach()

# preprocessor definitions
add_definitions(-D_USE_MATH_DEFINES -DGL_SILENCE_DEPRECATION)
if(WIN32)
	add_definitions(-DNOMINMAX)
endif()

foreach(lib ${VEGA_LIBRARIES})
	add_subdirectory(libraries/${lib})
endforeach()


# link everything together into one library
set(OBJECT_FILES)
foreach(lib ${VEGA_LIBRARIES})
	list(APPEND OBJECT_FILES $<TARGET_OBJECTS:${lib}>)
endforeach()

add_library(vega ${OBJECT_FILES})


# build utilities
if(VEGA_BUILD_UTILITIES)
	add_subdirectory(utilities)
endif()
