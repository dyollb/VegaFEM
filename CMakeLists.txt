cmake_minimum_required(VERSION 3.18)
project(VegaFEM CXX)
set(CMAKE_CXX_STANDARD 14)


# configuration options
option(VEGAFEM_USE_TBB "Build VegaFEM using TBB" ON)
option(VEGAFEM_USE_GLEW "Build OpenGL shaders using GLEW" OFF)
option(VEGAFEM_USE_GLUT "Build components that depend on GLUT" OFF)
option(VEGAFEM_USE_CG_TOOLKIT "Build GPU components using Cg Toolkit" OFF)
option(VEGAFEM_USE_CGAL "Build meshing components using CGAL" OFF)


# thirdparties
set(VEGAFEM_DEPENDENCIES)

find_package(OpenGL REQUIRED)
list(APPEND VEGAFEM_DEPENDENCIES OpenGL::GL)

if(VEGAFEM_USE_TBB)
	find_package(TBB REQUIRED)
	list(APPEND VEGAFEM_DEPENDENCIES ${TBB_IMPORTED_TARGETS})
	add_definitions(-DUSE_TBB)
endif()

if(VEGAFEM_USE_GLUT)
	find_package(GLUT REQUIRED)
	list(APPEND VEGAFEM_DEPENDENCIES GLUT::GLUT)
endif()

if(VEGAFEM_USE_GLEW)
	find_package(GLEW REQUIRED)
	list(APPEND VEGAFEM_DEPENDENCIES GLEW::GLEW) # glew_s for static
endif()

if(VEGAFEM_USE_CG_TOOLKIT)
	find_package(Cg QUIET) # todo: this does not work yet
endif()

if(VEGAFEM_USE_CGAL)
	find_package(CGAL REQUIRED)
	find_path(IGL_INCLUDE_DIRS igl/harmonic.h PATHS $ENV{IGL_ROOT}/include REQUIRED)
endif()

set(ENABLE_BLAS95 ON)
set(ENABLE_LAPACK95 ON)
find_package(MKL QUIET)
if(MKL_FOUND)
	set(VEGA_USE_INTEL_MKL ON)
	list(APPEND VEGAFEM_DEPENDENCIES MKL::MKL) # ${MKL_CORE} MKL::${MKL_BLAS95} MKL::${MKL_LAPACK95}
	get_target_property(MKL_INCLUDE_DIRS MKL::MKL INTERFACE_INCLUDE_DIRECTORIES)
	#include_directories(${MKL_INCLUDE_DIRS})
else()
	find_package(BLAS REQUIRED)
	list(APPEND VEGAFEM_DEPENDENCIES BLAS::BLAS)
	if(NOT APPLE)
		#find_path(BLAS_INCLUDE_DIRS cblas.h PATHS /usr/include /usr/local/include REQUIRED)
		#set_property(TARGET nana_lib PROPERTY INTERFACE_INCLUDE_DIRECTORIES
	endif()
	find_package(LAPACK REQUIRED)
	list(APPEND VEGAFEM_DEPENDENCIES LAPACK::LAPACK)
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

if(VEGAFEM_USE_GLUT)
	list(APPEND VEGA_LIBRARIES
		camera
		lighting
		openGLHelper
		renderVolumetricMesh
		sceneObject
		sceneObjectReduced
	)
endif()

if(VEGAFEM_USE_GLUT AND VEGAFEM_USE_GLEW)
	list(APPEND VEGA_LIBRARIES 
		glslPhong
	)
endif()

if(VEGAFEM_USE_CGAL)
	list(APPEND VEGA_LIBRARIES
		exactArithmetic
		libiglInterface
		virtualTets
	)
endif()

if(VEGAFEM_USE_CG_TOOLKIT)
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
target_link_libraries(vega ${VEGAFEM_DEPENDENCIES})

# build utilities
add_subdirectory(utilities)