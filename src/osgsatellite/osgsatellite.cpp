/*
*	This demo is a satellite image simulator
*/
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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


std::string readShaderFile(const std::string& filePath) {
    std::ifstream shaderFile(filePath);
    std::stringstream shaderStream;

    if (shaderFile.is_open()) {
        shaderStream << shaderFile.rdbuf();
        shaderFile.close();
        return shaderStream.str();
    }
    else {
        std::cerr << "Failed to open shader file: " << filePath << std::endl;
        return "";
    }
}

int main()
{
    osgViewer::Viewer* viewer = new osgViewer::Viewer;
    osg::ref_ptr<osg::Group> rootGroup = new osg::Group;
    osg::ref_ptr<osg::Group> SatelliteGroup = new osg::Group;// 前(x)右(y)下(z)，转序就是zyx(321)
    osg::ref_ptr<osg::Node> axesNode = osgDB::readNodeFile("axes.osgt");
    osg::ref_ptr<osg::MatrixTransform> mtScaleNode = new osg::MatrixTransform;
    mtScaleNode->setMatrix(osg::Matrix::scale(10, 10, 10));
    mtScaleNode->addChild(axesNode);

    osg::ref_ptr<osg::MatrixTransform> satellitePosMt = new osg::MatrixTransform;
    satellitePosMt->addChild(osgDB::readNodeFile("GPS.ive"));
    SatelliteGroup->addChild(satellitePosMt.get());


    osg::ref_ptr<osg::Geode> satellitePointGeode = new osg::Geode();
    osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
    geom->setDataVariance(osg::Object::DYNAMIC);
    geom->setUseVertexBufferObjects(true);
    // Define the vertex array
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->clear();
    vertices->push_back(osg::Vec3(0.0f, 0.0f, 0.0f)); // Point position
    vertices->push_back(osg::Vec3(5.0f, 0.0f, 0.0f)); // Point position
    vertices->push_back(osg::Vec3(10.0f, 0.0f, 0.0f)); // Point position
    geom->setVertexArray(vertices.get());
    osg::ref_ptr<osg::Vec4Array> colorArray = new osg::Vec4Array;
    colorArray->push_back(osg::Vec4(0.0f, 255.0f, 0.0f, 1.0f));
    colorArray->push_back(osg::Vec4(0.0f, 255.0f, 0.0f, 0.0f));
    geom->setColorArray(colorArray.get(), osg::Array::BIND_PER_VERTEX);

    // Set up the shader program
    osg::ref_ptr<osg::Program> program = new osg::Program;
    
    program->addShader(new osg::Shader(osg::Shader::VERTEX, readShaderFile("D:\\DevGitHub\\OSG-OE-TechDemos\\src\\osgsatellite\\satellitepoint.vs").c_str())); // Load the vertex shader
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, readShaderFile("D:\\DevGitHub\\OSG-OE-TechDemos\\src\\osgsatellite\\satellitepoint.fs").c_str())); // Load the fragment shader

    // Attach the shader program to the geometry
    geom->getOrCreateStateSet()->setAttributeAndModes(program.get());
    // Get the width of the window
    //unsigned int width = viewer->getCamera()->getViewport()->width();

    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_width_half", 1280)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_height_half", 720)); // for temporary
    geom->getOrCreateStateSet()->setMode(GL_POINT_SMOOTH, osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED);
    geom->getOrCreateStateSet()->setMode(GL_POINT_SPRITE, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
    geom->getOrCreateStateSet()->setAttributeAndModes(new osg::BlendFunc(), osg::StateAttribute::ON);
    // Set up point rendering
    geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, 2)); // don't know why, if count == 1, there will be nothing on the screen.
    // Enable point size attribute
    osg::ref_ptr<osg::Point> pointAttrib = new osg::Point;
    pointAttrib->setSize(12.0); // Set the size of the point
    geom->getOrCreateStateSet()->setAttributeAndModes(pointAttrib.get());

    satellitePointGeode->addDrawable(geom.get());
    satellitePointGeode->getOrCreateStateSet()->setRenderBinDetails(1, "RenderBin");
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Depth(), osg::StateAttribute::OFF);
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Light(), osg::StateAttribute::OFF);
    satellitePosMt->addChild(satellitePointGeode.get());
    rootGroup->addChild(SatelliteGroup.get());
    viewer->setSceneData(rootGroup.get());
    return viewer->run();
}

