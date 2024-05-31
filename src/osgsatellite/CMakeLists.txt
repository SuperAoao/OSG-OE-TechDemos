add_executable(osgSatellite osgsatellite.cpp satellitepoint.vs satellitepoint.fs)

target_link_libraries(osgSatellite ${OPENSCENEGRAPH_LIBRARIES})