include(${CMAKE_SOURCE_DIR}/cmake/CopyFilesUtility.cmake)

set(MY_TARGET_NAME osgSatellite)

set(COMMON_SRC
${PROJECT_SOURCE_DIR}/src/public_header/helper.h
${PROJECT_SOURCE_DIR}/src/public_source/helper.cpp
${PROJECT_SOURCE_DIR}/src/public_header/pbrstateset.h
${PROJECT_SOURCE_DIR}/src/public_source/pbrstateset.cpp
)

set(COMMON_SHADER_SRC
${PROJECT_SOURCE_DIR}/src/public_shader/pbr.vs
${PROJECT_SOURCE_DIR}/src/public_shader/pbr.gs
${PROJECT_SOURCE_DIR}/src/public_shader/pbr.fs
)

set(SHADER_SRC
satellitepoint.vs
satellitepoint.fs
)

source_group("Common Files" FILES ${COMMON_SRC})

source_group("Common Shader Files" FILES ${COMMON_SHADER_SRC})

source_group("Shader Files" FILES ${SHADER_SRC})


include_directories(${PROJECT_SOURCE_DIR}/src/public_header)

add_executable(${MY_TARGET_NAME} osgsatellite.cpp ${COMMON_SRC} ${COMMON_SHADER_SRC} ${SHADER_SRC})

target_link_libraries(${MY_TARGET_NAME} ${OPENSCENEGRAPH_LIBRARIES})

# copy shader files

set(SHADER_DEST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders")

copy_files_on_build(${MY_TARGET_NAME} ${SHADER_DEST_DIR}
${COMMON_SHADER_SRC}
${SHADER_SRC}
)

#generator expression $<CONFIG>

set_target_properties(${MY_TARGET_NAME} PROPERTIES 
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>"
)