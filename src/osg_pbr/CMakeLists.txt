include(${CMAKE_SOURCE_DIR}/cmake/CopyFilesUtility.cmake)

set(MY_TARGET_NAME osgPBR)

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

source_group("Common Files" FILES ${COMMON_SRC})

source_group("Common Shader" FILES ${COMMON_SHADER_SRC})


add_executable(osgPBR main.cpp ${COMMON_SRC} ${COMMON_SHADER_SRC})

include_directories(${PROJECT_SOURCE_DIR}/src/public_header)

target_link_libraries(${MY_TARGET_NAME} ${OPENSCENEGRAPH_LIBRARIES})

# copy shader files

set(SHADER_DEST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders")

copy_files_on_build(${MY_TARGET_NAME} ${SHADER_DEST_DIR}
${COMMON_SHADER_SRC}
)

#generator expression $<CONFIG>

set_target_properties(${MY_TARGET_NAME} PROPERTIES 
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>"
)
