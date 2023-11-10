if(SPOOLES_INCLUDE_DIR AND SPOOLES_LIBRARY)
  set(SPOOLES_FIND_QUIETLY TRUE)
endif()

if(BUILD_SHARED_LIBS)
  set(_spooles_lib_name spooles.so)
else()
  set(_spooles_lib_name spooles.a)
endif()

find_path(SPOOLES_INCLUDE_DIR LinSol/Bridge.h
          PATHS ${SPOOLES_ROOT}/include /usr/include /usr/local/include
                REQUIRED)
find_library(
  SPOOLES_LIBRARY
  NAMES spooles ${_spooles_lib_name}
  PATHS ${SPOOLES_ROOT}/lib /usr/lib /usr/local/lib REQUIRED)

add_library(SPOOLES INTERFACE)
target_include_directories(SPOOLES INTERFACE ${SPOOLES_INCLUDE_DIR})
target_link_libraries(SPOOLES INTERFACE ${SPOOLES_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if all
# listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPOOLES DEFAULT_MSG SPOOLES_INCLUDE_DIR
                                  SPOOLES_LIBRARY)

mark_as_advanced(SPOOLES_INCLUDE_DIR SPOOLES_LIBRARY)
