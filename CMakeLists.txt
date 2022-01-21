cmake_minimum_required(VERSION 3.1)
project(VegaFEM CXX)
set(CMAKE_CXX_STANDARD 14)

# find dependencies
find_package(OpenGL REQUIRED)

find_package(GLUT QUIET)

find_package(GLEW QUIET)

find_package(Cg QUIET)

set(BLA_VENDOR Intel10_64lp)
find_package(BLAS QUIET)
if(BLAS_FOUND)
	set(USE_INTEL_MKL ON)
else()
	set(USE_INTEL_MKL OFF)
    unset(BLA_VENDOR)
    find_package(BLAS REQUIRED)
	find_path(BLAS_INCLUDE_DIRS cblas.h PATHS
		/usr/include
		/usr/local/include
		$ENV{BLAS_HOME}/include REQUIRED)
endif()
if(BLAS_FOUND)
	find_package(LAPACK REQUIRED)
endif(BLAS_FOUND)

find_package(CGAL QUIET)

find_package(TBB QUIET)

find_path(IGL_INCLUDE_DIRS igl/harmonic.h)


include_directories(
	${OPENGL_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include
	${CMAKE_BINARY_DIR}
)

if (GLUT_FOUND)
	include_directories(${GLUT_INCLUDE_DIR})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libraries/include/vega-config.h.in ${CMAKE_BINARY_DIR}/vega-config.h)


if (TBB_FOUND)
	get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
	add_definitions(-DUSE_TBB)
	include_directories(${TBB_INCLUDE_DIRS})
endif()

if (BLAS_INCLUDE_DIRS)
	include_directories(${BLAS_INCLUDE_DIRS})
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
		sparseSolver
	)
	if(OPENGL_FOUND AND GLUT_FOUND)
		list(APPEND VEGA_CORE_LIBRARIES sceneObjectReduced)
	endif()
endif()

if(CGAL_FOUND AND IGL_INCLUDE_DIRS)
	list(APPEND VEGA_CORE_LIBRARIES
		exactArithmetic
		libiglInterface
		virtualTets
	)
endif()


foreach(lib ${VEGA_CORE_LIBRARIES})
	include_directories(libraries/${lib})
endforeach()

add_definitions(-D_USE_MATH_DEFINES -DGL_SILENCE_DEPRECATION)
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