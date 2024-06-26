# Locate zlib
# This module defines
# ZLIB_LIBRARY
# ZLIB_FOUND, if false, do not try to link to zlib 
# ZLIB_INCLUDE_DIR, where to find the headers
#
# $ZLIB_DIR is an environment variable that would
# correspond to the ./configure --prefix=$ZLIB_DIR
# used in building zlib.
#
# Created by Ulrich Hertlein. 

# prefer FindZLIB from cmake distribution

#HGT 2020.11.10 x00038 使用cmake自带的findzlib 会导致无法打开输入文件optimized.lib报错
#if(EXISTS ${CMAKE_ROOT}/Modules/FindZLIB.cmake)
#  include(${CMAKE_ROOT}/Modules/FindZLIB.cmake)
#
#  if(ZLIB_FOUND)
#    return()
#  endif()
#endif()
#end


FIND_PATH(ZLIB_INCLUDE_DIR zlib.h
    $ENV{ZLIB_DIR}/include
    $ENV{ZLIB_DIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/include
    /usr/include
    /sw/include # Fink
    /opt/local/include # DarwinPorts
    /opt/csw/include # Blastwave
    /opt/include
    /usr/freeware/include
)

FIND_LIBRARY(ZLIB_LIBRARY 
    NAMES z libz zlib
    PATHS
    $ENV{ZLIB_DIR}/lib
    $ENV{ZLIB_DIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/lib
    /usr/lib
    /sw/lib
    /opt/local/lib
    /opt/csw/lib
    /opt/lib
    /usr/freeware/lib64
)

SET(ZLIB_FOUND "NO")
IF(ZLIB_LIBRARY AND ZLIB_INCLUDE_DIR)
    SET(ZLIB_FOUND "YES")
ENDIF(ZLIB_LIBRARY AND ZLIB_INCLUDE_DIR)


