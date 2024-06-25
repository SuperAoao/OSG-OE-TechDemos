/*
* This demo illustrates a physically-based renderer implemented in OSG
*/
#include "helper.h"
#include "pbrstateset.h"

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

// We update the uniform variable osg_ViewMatrixInverse in a customized OSG SceneView, so we won’t use this callback in this demo.
class CameraUpdateCallback : public osg::NodeCallback {
public:
    CameraUpdateCallback(osg::StateSet* stateSet) : m_stateSet(stateSet) {};
    virtual void operator()(osg::Node* node, osg::NodeVisitor* nv) {
        osg::Camera* camera = dynamic_cast<osg::Camera*>(node);
        if (camera) {
            osg::Matrix viewMatrix = camera->getViewMatrix();
            osg::Matrix inverseViewMatrix = osg::Matrix::inverse(viewMatrix);

            osg::Uniform* uniform = m_stateSet->getOrCreateUniform("osg_ViewMatrixInverse", osg::Uniform::FLOAT_MAT4);
            uniform->set(inverseViewMatrix);
        }

        traverse(node, nv);
    }
private:
    osg::StateSet* m_stateSet;
};

// Update light related variables 
class PBRNodeUpdateCallback : public osg::NodeCallback {
public:
    PBRNodeUpdateCallback(osgViewer::Viewer* viewer, osg::MatrixTransform* mtTrans, osg::StateSet* pbrStateSet) : m_viewer(viewer),
        m_mtTrans(mtTrans),
        m_pbrStateSet(pbrStateSet)
    {};
    virtual void operator()(osg::Node* node, osg::NodeVisitor* nv) 
    {
        osg::Node* pbrNode = node;
        osg::Vec4 ltDir(1, 0, 0, 0);
        osg::Vec3 ltCol(60, 60, 60);
        if (m_viewer)
        {
            osg::Light* lt = m_viewer->getLight();
            ltDir = lt->getPosition();
        }

        osg::Uniform* lightNumUni = m_pbrStateSet->getOrCreateUniform("globalLightNum", osg::Uniform::INT);
        if (lightNumUni)
        {
            int lightNum = 0;
            lightNumUni->get(lightNum);
            if (lightNum == 0)
            {
                lightNumUni->set(1);
            }
        }
  
        osg::Uniform* lightPosUni = m_pbrStateSet->getOrCreateUniform("globalLightPositions", osg::Uniform::FLOAT_VEC4, 200.0f);
        if (lightPosUni)
        {
            lightPosUni->setElement(0, ltDir);
        }
        
        osg::Uniform* lightColUni = m_pbrStateSet->getOrCreateUniform("globalLightColors", osg::Uniform::FLOAT_VEC3, 200.0f);
        if (lightColUni)
        {
            lightColUni->setElement(0, ltCol);
        }
        // Continue traversal to children
        traverse(node, nv);
    }
private:
    osgViewer::Viewer* m_viewer;
    osg::MatrixTransform* m_mtTrans;
    osg::StateSet* m_pbrStateSet;
};

int main()
{
    osgViewer::Viewer* viewer = new osgViewer::Viewer;
    osg::ref_ptr<osg::Group> rootGroup = new osg::Group;
    osg::ref_ptr<osg::Group> SatelliteGroup = new osg::Group;
    osg::ref_ptr<osg::Node> axesNode = osgDB::readNodeFile("axes.osgt");//createCoordinateAxes();
    axesNode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF | osg::StateAttribute::PROTECTED);
    osg::ref_ptr<osg::MatrixTransform> axesScaleMt = new osg::MatrixTransform;
    axesScaleMt->setMatrix(osg::Matrix::scale(10, 10, 10));
    axesScaleMt->addChild(axesNode);

    osg::ref_ptr<osg::MatrixTransform> satellitePosMt = new osg::MatrixTransform;
    //satellitePosMt->addChild(osgDB::readNodeFile("GPS.ive"));
    //osg::Node* modelNode = osgDB::readNodeFile("D:\\Dev\\ht5\\install\\FreeXGISDesktop-Data-FXG_3.8.2\\Data\\Model\\ive\\weixing.ive");
    osg::Node* modelNode = osgDB::readNodeFile("E:\\3DAssets\\shield\\Low.fbx");
    satellitePosMt->addChild(modelNode);
    // PBR
    osg::ref_ptr<PBRStateSet> pbrStateSet = new PBRStateSet();
    pbrStateSet->loadTextures("E:\\3DAssets\\shield\\textur", PBRStateSet::PNG);
    //PBRStateSet::setEnvironmentLighting(osg::Vec3f(-10.0f, 0.0f, 10.0f));
    modelNode->setStateSet(pbrStateSet.get());
    modelNode->setUpdateCallback(new PBRNodeUpdateCallback(viewer, satellitePosMt.get(), pbrStateSet.get()));
    SatelliteGroup->addChild(satellitePosMt.get());
   
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
    //camera->setUpdateCallback(new CameraUpdateCallback(pbrStateSet.get()));
    
    viewer->getLight()->setPosition(osg::Vec4(1.0, 0.0, 1.0, 0.0));
    viewer->getLight()->setAmbient(osg::Vec4(1.0, 1.0, 1.0, 1.0));

    return viewer->run();
}

