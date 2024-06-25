cmake_minimum_required(VERSION 3.27)

project(OSG-OE-TechDemos)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(CMAKE_DEBUG_POSTFIX "d" CACHE STRING "add a postfix, usually d on windows")
set(CMAKE_RELEASE_POSTFIX "" CACHE STRING "add a postfix, usually empty on windows")
set(CMAKE_RELWITHDEBINFO_POSTFIX "" CACHE STRING "add a postfix, usually empty on windows")
set(CMAKE_MINSIZEREL_POSTFIX "" CACHE STRING "add a postfix, usually empty on windows")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")

find_package(OpenGL)
find_package(JPEG)
find_package(PNG)
find_package(CURL)
find_package(ZLIB)
find_package(OpenSceneGraph 3.4.0 REQUIRED 
	osg 
	osgGA 
	osgDB 
	osgViewer 
	osgUtil 
	osgText 
	osgQt 
	osgVolume 
	osgTerrain 
	osgSim 
	osgFX 
	osgAnimation 
	osgManipulator
	osgParticle)
message("OPENSCENEGRAPH_FOUND: ${OPENSCENEGRAPH_FOUND}" )
message("OPENSCENEGRAPH_VERSION: ${OPENSCENEGRAPH_VERSION}")
message("OPENSCENEGRAPH_INCLUDE_DIRS: ${OPENSCENEGRAPH_INCLUDE_DIRS}")

include_directories(${OPENSCENEGRAPH_INCLUDE_DIRS})

# src
add_subdirectory(src/osg_satellite)
add_subdirectory(src/osg_pbr)


