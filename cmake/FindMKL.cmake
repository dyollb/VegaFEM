# ##############################################################################
#
# \file      cmake/FindMKL.cmake \author    J. Bakosi \copyright 2012-2015,
# Jozsef Bakosi, 2016, Los Alamos National Security, LLC. \brief     Find the
# Math Kernel Library from Intel \date      Thu 26 Jan 2017 02:05:50 PM MST
#
# ##############################################################################

# Find the Math Kernel Library from Intel
#
# MKL_FOUND - System has MKL MKL_INCLUDE_DIRS - MKL include files directories
# MKL_LIBRARIES - The MKL libraries MKL_INTERFACE_LIBRARY - MKL interface
# library MKL_SEQUENTIAL_LAYER_LIBRARY - MKL sequential layer library
# MKL_CORE_LIBRARY - MKL core library
#
# The environment variables MKLROOT and INTEL are used to find the library.
# Everything else is ignored. If MKL is found "-DMKL_LP64" is added to
# CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
#
# Example usage:
#
# find_package(MKL) if(MKL_FOUND) target_link_libraries(TARGET ${MKL_LIBRARIES})
# endif()

# If already in cache, be silent
if(MKL_INCLUDE_DIRS
   AND MKL_LIBRARIES
   AND MKL_INTERFACE_LIBRARY
   AND MKL_SEQUENTIAL_LAYER_LIBRARY
   AND MKL_RT_LIBRARY
   AND MKL_CORE_LIBRARY)
  set(MKL_FIND_QUIETLY TRUE)
endif()

if(NOT BUILD_SHARED_LIBS)
  if(WIN32)
    set(INT_LIB "mkl_intel_lp64.lib")
    set(SEQ_LIB "mkl_sequential.lib")
    set(THR_LIB "mkl_intel_thread.lib")
    set(COR_LIB "mkl_core.lib")
    set(RT_LIB "mkl_rt.lib")
  else()
    set(INT_LIB "libmkl_intel_lp64.a")
    set(SEQ_LIB "libmkl_sequential.a")
    set(THR_LIB "libmkl_intel_thread.a")
    set(COR_LIB "libmkl_core.a")
    set(RT_LIB "libmkl_rt.a")
  endif()
else()
  set(INT_LIB "mkl_intel_lp64")
  set(SEQ_LIB "mkl_sequential")
  set(THR_LIB "mkl_intel_thread")
  set(COR_LIB "mkl_core")
  set(RT_LIB "mkl_rt")
endif()

if(MSVC)
  set(ProgramFilesx86 "ProgramFiles(x86)")
  set(INTEL_ROOT_DEFAULT
      $ENV{${ProgramFilesx86}}/IntelSWTools/compilers_and_libraries/windows)
else()
  set(INTEL_ROOT_DEFAULT "/opt/intel")
endif()
set(INTEL_ROOT
    ${INTEL_ROOT_DEFAULT}
    CACHE PATH "Folder contains intel libs")
find_path(
  MKL_ROOT include/mkl.h
  PATHS $ENV{MKL_ROOT} $ENV{MKLROOT} ${INTEL_ROOT}/mkl
  DOC "Folder contains MKL")

find_path(
  MKL_INCLUDE_DIR
  NAMES mkl.h
  HINTS ${MKL_ROOT}/include /usr/include/mkl)

find_library(
  MKL_INTERFACE_LIBRARY
  NAMES ${INT_LIB}
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64 ${MKL_ROOT}/lib/intel64_win
        ${INTEL_ROOT}/mkl/lib/intel64)

find_library(
  MKL_SEQUENTIAL_LAYER_LIBRARY
  NAMES ${SEQ_LIB}
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64 ${INTEL_ROOT}/mkl/lib/intel64)

find_library(
  MKL_THREADED_LAYER_LIBRARY
  NAMES ${THR_LIB}
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64 ${INTEL_ROOT}/mkl/lib/intel64)

find_library(
  MKL_CORE_LIBRARY
  NAMES ${COR_LIB}
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64 ${INTEL_ROOT}/mkl/lib/intel64)
  
find_library(
  MKL_RT_LIBRARY
  NAMES ${RT_LIB}
  PATHS ${MKL_ROOT}/lib ${MKL_ROOT}/lib/intel64 ${INTEL_ROOT}/mkl/lib/intel64)
      
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
if (MKL_RUNTIME_LAYERS)
    message(STATUS "linking against ${MKL_RT_LIBRARY}")
    set(MKL_LIBRARIES ${MKL_RT_LIBRARY})
else()
    set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADED_LAYER_LIBRARY}
                  ${MKL_CORE_LIBRARY})
endif()

if(NOT WIN32 AND NOT APPLE)
  # Added -Wl block to avoid circular dependencies.
  # https://stackoverflow.com/questions/5651869/what-are-the-start-group-and-end-group-command-line-options
  # https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
  set(MKL_LIBRARIES -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group)
elseif(APPLE)
  # MacOS does not support --start-group and --end-group
  set(MKL_LIBRARIES -Wl,${MKL_LIBRARIES} -Wl,)
endif()

# message("1 ${MKL_INCLUDE_DIR}") message("2 ${MKL_INTERFACE_LIBRARY}")
# message("3 ${MKL_SEQUENTIAL_LAYER_LIBRARY}") message("4 ${MKL_CORE_LIBRARY}")

if(MKL_INCLUDE_DIR
   AND MKL_INTERFACE_LIBRARY
   AND MKL_SEQUENTIAL_LAYER_LIBRARY
   AND MKL_THREADED_LAYER_LIBRARY
   AND MKL_RT_LIBRARY   
   AND MKL_CORE_LIBRARY)

  if(NOT DEFINED ENV{CRAY_PRGENVPGI}
     AND NOT DEFINED ENV{CRAY_PRGENVGNU}
     AND NOT DEFINED ENV{CRAY_PRGENVCRAY}
     AND NOT DEFINED ENV{CRAY_PRGENVINTEL}
     AND NOT MSVC)
    set(ABI "-m64")
  endif()

  # set(MKL_DEFINITIONS "-DMKL_LP64 ${ABI}")

  add_library(MKL INTERFACE)
  add_library(MKL::MKL ALIAS MKL)
  target_include_directories(MKL INTERFACE ${MKL_INCLUDE_DIR})
  target_compile_definitions(MKL INTERFACE ${MKL_DEFINITIONS})
  target_link_libraries(MKL INTERFACE ${MKL_LIBRARIES})

else()
  set(MKL_INCLUDE_DIRS "")
  set(MKL_LIBRARIES "")
  set(MKL_INTERFACE_LIBRARY "")
  set(MKL_SEQUENTIAL_LAYER_LIBRARY "")
  set(MKL_THREADED_LAYER_LIBRARY "")
  set(MKL_RT_LIBRARY "")
  set(MKL_CORE_LIBRARY "")

endif()

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if all
# listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MKL
  DEFAULT_MSG
  MKL_LIBRARIES
  MKL_INCLUDE_DIRS
  MKL_INTERFACE_LIBRARY
  MKL_THREADED_LAYER_LIBRARY
  MKL_RT_LIBRARY
  MKL_CORE_LIBRARY)

mark_as_advanced(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY
                 MKL_THREADED_LAYER_LIBRARY MKL_CORE_LIBRARY MKL_RT_LIBRARY)
