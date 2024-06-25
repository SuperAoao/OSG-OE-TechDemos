# CopyFilesUtility.cmake

if(NOT COMMAND copy_files_on_build)
    function(copy_files_on_build TARGET_NAME DEST_DIR)
        # Get the list of files to copy (all arguments after DEST_DIR)
        set(FILES_TO_COPY ${ARGN})

        # Create the destination directory if it doesn't exist
        if(NOT EXISTS ${dir})
            file(MAKE_DIRECTORY ${DEST_DIR})
        endif()
        

        # Create a custom command to copy each file
        foreach(FILE_TO_COPY ${FILES_TO_COPY})
            get_filename_component(FILE_NAME ${FILE_TO_COPY} NAME)
            add_custom_command(
                TARGET ${TARGET_NAME}
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                ${FILE_TO_COPY}
                ${DEST_DIR}/${FILE_NAME}
                COMMENT "Copying ${FILE_TO_COPY} to ${DEST_DIR}"
            )
        endforeach()
    endfunction()
endif()