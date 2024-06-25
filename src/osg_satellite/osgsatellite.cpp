/*
*	This demo is a satellite image simulator
*/
#include "helper.h"
#include "starimgsimutils.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include <osgViewer/Viewer>
#include <osg/Group>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/vec3>
#include <osg/Vec3f>
#include <osgGA/StateSetManipulator>
#include <osg/MatrixTransform>
#include <osgViewer/ViewerEventHandlers>
#include <osg/StateSet>
#include <osg/BlendFunc>
#include <osg/Depth>
#include <osg/Point>
#include <osg/Uniform>
#include <osgDB/Registry>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osg/Image>
#include <osg/NodeVisitor>
#include <osg/Drawable>
#include <osg/LineWidth>
#include <osg/Light>
#include <osg/LightSource>
#include <osg/LOD>
#include <osgViewer/ViewerEventHandlers>

osg::ref_ptr<osg::Geode> createCoordinateAxes()
{
    osg::ref_ptr<osg::Geode> geode = new osg::Geode();

    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();

    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->push_back(osg::Vec3(0, 0, 0));
    vertices->push_back(osg::Vec3(1, 0, 0));
    vertices->push_back(osg::Vec3(0, 0, 0));
    vertices->push_back(osg::Vec3(0, 1, 0));
    vertices->push_back(osg::Vec3(0, 0, 0));
    vertices->push_back(osg::Vec3(0, 0, 1));

    geometry->setVertexArray(vertices);

    osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(osg::PrimitiveSet::LINES, 0);
    indices->push_back(0);
    indices->push_back(1);
    indices->push_back(2);
    indices->push_back(3);
    indices->push_back(4);
    indices->push_back(5);

    geometry->addPrimitiveSet(indices);

    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0)); // Red for X-axis
    colors->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0)); // Red for X-axis
    colors->push_back(osg::Vec4(0.0, 1.0, 0.0, 1.0)); // Green for Y-axis
    colors->push_back(osg::Vec4(0.0, 1.0, 0.0, 1.0)); // Green for Y-axis
    colors->push_back(osg::Vec4(0.0, 0.0, 1.0, 1.0)); // Blue for Z-axis
    colors->push_back(osg::Vec4(0.0, 0.0, 1.0, 1.0)); // Blue for Z-axis

    geometry->setColorArray(colors, osg::Array::BIND_PER_VERTEX);

    osg::ref_ptr<osg::LineWidth> linewidth = new osg::LineWidth();
    linewidth->setWidth(5.0f);

    geometry->getOrCreateStateSet()->setAttributeAndModes(linewidth, osg::StateAttribute::ON);

    geode->addDrawable(geometry);

    return geode;
}

int main()
{
    // calculate irradiance of satellite 
    //double E = computeIrradianceOfSatellite(648610.0, 0.25, 0.290);
    double E = StarImgSimUtils::computeIrradianceOfSatellite(20440.0, StarImgSimUtils::g_rho, 0.800);
    double target_mag = StarImgSimUtils::irradianceToApparentMagnitude(E);
    double d_len = 0.03; // unit m
    double t_exposure = 5;    // unit s
    double n_pe = StarImgSimUtils::getPhotonsFromMagnitude(target_mag, d_len, t_exposure);

    osg::ref_ptr<osg::Texture2D> PSFTex = StarImgSimUtils::genPSFTexture(12, 0.85);
    // apply energy
    double grayRatio = n_pe / StarImgSimUtils::g_n_e;
    unsigned char* imgData = PSFTex->getImage()->data();
    for (int i = 0; i < 12; ++i)
    {
        for (int j = 0; j < 12; ++j)
        {
            double proportion = (double)imgData[i * 12 + j] / 256.0;
            double grayVal = std::min(proportion* grayRatio,255.0);
            imgData[i * 12 + j] = (unsigned char)(grayVal);
        }
    }

    osgViewer::Viewer* viewer = new osgViewer::Viewer;
    osg::ref_ptr<osg::Group> rootGroup = new osg::Group;
    osg::ref_ptr<osg::Group> SatelliteGroup = new osg::Group;// 前(x)右(y)下(z)，转序就是zyx(321)
    osg::ref_ptr<osg::Node> axesNode = osgDB::readNodeFile("axes.osgt");//createCoordinateAxes();
    axesNode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED);
    osg::ref_ptr<osg::MatrixTransform> axesScaleMt = new osg::MatrixTransform;
    axesScaleMt->setMatrix(osg::Matrix::scale(10, 10, 10));
    axesScaleMt->addChild(axesNode);

    osg::ref_ptr<osg::MatrixTransform> satellitePosMt = new osg::MatrixTransform;
    // Non PBR
    osg::ref_ptr<osg::Group> satelliteNonPBRGroup = new osg::Group;
    osg::ref_ptr<osg::MatrixTransform> satelliteNonPBRMT = new osg::MatrixTransform;
    satelliteNonPBRMT->setMatrix(osg::Matrix::rotate(osg::DegreesToRadians(90.0), osg::Vec3d(0, 0, 1),
        osg::DegreesToRadians(-180.0), osg::Vec3d(1, 0, 0),
        osg::DegreesToRadians(0.0), osg::Vec3d(0, 0, 0)));

    satelliteNonPBRGroup->addChild(satelliteNonPBRMT.get());
    satelliteNonPBRMT->addChild(osgDB::readNodeFile("GPS.ive"));
    //satelliteNonPBRMT->addChild(osgDB::readNodeFile("Bd1h.ive"));



    satellitePosMt->addChild(satelliteNonPBRGroup.get());
    SatelliteGroup->addChild(satellitePosMt.get());
    SatelliteGroup->addChild(axesScaleMt.get()); // add axesScaleMt to scene graph

    osg::ref_ptr<osg::Geode> satellitePointGeode = new osg::Geode();
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
    geom->setDataVariance(osg::Object::DYNAMIC);
    geom->setUseVertexBufferObjects(true);
    // Define the vertex array
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->clear();
    vertices->push_back(osg::Vec3(0.0f, 0.0f, 0.0f)); // Point position
    //vertices->push_back(osg::Vec3(500000.0f, 0.0f, 0.0f)); // Point position
    //vertices->push_back(osg::Vec3(10.0f, 0.0f, 0.0f)); // Point position
    geom->setVertexArray(vertices.get());
    geom->dirtyBound();
    /*osg::ref_ptr<osg::Vec4Array> colorArray = new osg::Vec4Array;
    colorArray->push_back(osg::Vec4(0.0f, 255.0f, 0.0f, 1.0f));
    colorArray->push_back(osg::Vec4(0.0f, 255.0f, 0.0f, 0.0f));
    geom->setColorArray(colorArray.get(), osg::Array::BIND_PER_VERTEX);*/

    // Set up the shader program
    osg::ref_ptr<osg::Program> program = new osg::Program;
    
    program->addShader(new osg::Shader(osg::Shader::VERTEX, Helper::readShaderFile("..\\shaders\\satellitepoint.vs").c_str())); // Load the vertex shader
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, Helper::readShaderFile("..\\shaders\\satellitepoint.fs").c_str())); // Load the fragment shader

    // Attach the shader program to the geometry
    geom->getOrCreateStateSet()->setAttributeAndModes(program.get());
    // Get the width of the window
    //unsigned int width = viewer->getCamera()->getViewport()->width();
    geom->getOrCreateStateSet()->setTextureAttributeAndModes(0, PSFTex.get(), osg::StateAttribute::ON);
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_width_half", 1280)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_height_half", 720)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("n_pe", n_pe));
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("n_e", StarImgSimUtils::g_n_e));
    geom->getOrCreateStateSet()->setMode(GL_POINT_SMOOTH, osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED);
    geom->getOrCreateStateSet()->setMode(GL_POINT_SPRITE, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
    geom->getOrCreateStateSet()->setAttributeAndModes(new osg::BlendFunc(), osg::StateAttribute::ON);
    // Set up point rendering
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, 1)); // don't know why, if count == 1, there will be nothing on the screen.
    // still don't know why, but finally it works out when I add updateDarkPoint func
    // Enable point size attribute
    osg::ref_ptr<osg::Point> pointAttrib = new osg::Point;
    pointAttrib->setSize(12.0); // Set the size of the point
    geom->getOrCreateStateSet()->setAttributeAndModes(pointAttrib.get()); 
    
    viewer->getLight()->setPosition(osg::Vec4(0.0, 0.0, 1.0,0.0));
    viewer->getLight()->setDirection(osg::Vec3(0.0, 0.0, -1.0));
    geom->setUpdateCallback(new StarImgSimUtils::DistanceUpdateCallback(viewer->getCamera(), viewer->getLight(),geom.get(), PSFTex.get()));
    satellitePointGeode->addDrawable(geom.get());
    satellitePointGeode->getOrCreateStateSet()->setRenderBinDetails(1, "RenderBin");
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Depth(), osg::StateAttribute::OFF);
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Light(), osg::StateAttribute::OFF);
    osg::ref_ptr<osg::LOD> pointLOD = new osg::LOD;
    pointLOD->addChild(satellitePointGeode.get(), 3000, 100000000);
    //satellitePosMt->addChild(satellitePointGeode.get());
    satellitePosMt->addChild(pointLOD.get());

    // Sun
    // Create a sphere to represent the sun
    osg::ref_ptr<osg::Sphere> sunShape = new osg::Sphere(osg::Vec3(0, 0, 0), 30.0f); // Center at origin, radius 1.0
    osg::ref_ptr<osg::ShapeDrawable> sunDrawable = new osg::ShapeDrawable(sunShape);

    // Optionally, set the color of the sun
    sunDrawable->setColor(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f)); // Bright yellow color

    // Create a geode to hold the drawable
    osg::ref_ptr<osg::Geode> sunGeode = new osg::Geode();
    sunGeode->addDrawable(sunDrawable.get());

    osg::ref_ptr<osg::MatrixTransform> sunMt = new osg::MatrixTransform;
    sunMt->setMatrix(osg::Matrix::translate(0.0f, 0.0f, 100.0f));
    sunMt->addChild(sunGeode.get());
    
    //rootGroup->addChild(sunMt.get()); // add Sun to scene graph
    rootGroup->addChild(SatelliteGroup.get());
    viewer->setSceneData(rootGroup.get());
    viewer->addEventHandler(new osgViewer::StatsHandler());
    double fovy, aspectRatio, zNear, zFar;
    viewer->getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, zFar);
    viewer->getCamera()->setCullingMode(osg::CullSettings::NO_CULLING);
    viewer->getCamera()->setProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, 100000);

    // Define the traits for the small window.
    osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits();
    traits->x = 100;           // Window position x
    traits->y = 100;           // Window position y
    traits->width = 800;       // Window width
    traits->height = 600;      // Window height
    traits->windowDecoration = true; // Enable window decorations (title bar, borders)
    traits->doubleBuffer = true;  // Enable double buffering

    // Create a graphics context based on the defined traits.
    osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

    if (!gc)
    {
        std::cerr << "Failed to create graphics context." << std::endl;
        return 1;
    }

    // Set up the camera for the viewer.
    osg::Camera* camera = viewer->getCamera();
    camera->setGraphicsContext(gc.get());
    camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));

    // Set the clear color for the window.
    camera->setClearColor(osg::Vec4(0.2f, 0.2f, 0.2f, 1.0f));

    // Set up projection matrix.
    camera->setProjectionMatrixAsPerspective(
        30.0f, static_cast<double>(traits->width) / static_cast<double>(traits->height), 1.0f, 100000.0f);

    return viewer->run();
}

