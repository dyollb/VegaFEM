macro(install_app exe)
set(options JUNK1)
set(oneValueArgs JUNK2)
set(multiValueArgs PLUGIN_NAMES LIBRARY_DIRS)
cmake_parse_arguments(INSTALL_APP "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

set(APP_NAME ${exe})
set(PLUGIN_NAMES ${INSTALL_APP_PLUGIN_NAMES})
set(BUNDLE_LIB_DIRS ${INSTALL_APP_LIBRARY_DIRS})
set(BUNDLE_DIR "${APP_NAME}")
set(LIBRARY_DIR ${BUNDLE_DIR})
if(APPLE)
    set(LIBRARY_DIR "${BUNDLE_DIR}/${APP_NAME}.app/Contents/Frameworks")
endif()

install(TARGETS ${APP_NAME} ${PLUGIN_NAMES}
    BUNDLE  DESTINATION ${BUNDLE_DIR} COMPONENT Runtime
    RUNTIME DESTINATION ${LIBRARY_DIR} COMPONENT Runtime
    LIBRARY DESTINATION ${LIBRARY_DIR} COMPONENT Runtime
)

set(BUNDLE_PLUGINS "")
foreach(plugin ${PLUGIN_NAMES})
    list(APPEND BUNDLE_PLUGINS 
        ${CMAKE_INSTALL_PREFIX}/${LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${plugin}${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
endforeach()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Bundle.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/Bundle.cmake
    @ONLY
)

install(SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/Bundle.cmake)

endmacro()
