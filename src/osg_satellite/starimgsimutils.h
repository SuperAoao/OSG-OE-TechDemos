#pragma once
#include <osg/vec3>
#include <osg/Texture2D>
#include <osg/Light>
#include <osg/Drawable>
#include <osg/Group>
#include <osg/NodeCallback>
#include <osgGA/GUIEventHandler>
#include <osgViewer/Viewer>
#include <osg/Material>
#include <osg/StateSet>

namespace StarImgSimUtils
{
// exposure time 1s, d_length 1 mm^2
const double g_CCDphotons = 19100;
// satellite
const double g_E_sun = 573; // Irradiance of the sun measured in[W/m^2]
const double g_sun_app_mag = -26.7; // sun apparent magnitude
const double g_rho = 0.25; // effective reflectance
// CCD gray parameters
const double g_n_ec = 45000;    // CCD per pixel max photons
const double g_n_e = g_n_ec / 256;    // photons per gray scale

// we consider center point (xc,yc) = (0,0)
double simplifiedGaussianFunc(const double x, const double y, const double xc, const double yc, const double sigma);

double getEnergyProportion(const double x_start, const double y_start, const double x_end, const double y_end, const double steps, const double sigma);

osg::Texture2D* genPSFTexture(unsigned int pixelSize, const double sigma);

// get photons from star magitude on CCD
double getPhotonsFromMagnitude(double mv, double d_len, double t_in);

/**
 * @brief
 * @param r -> the distance between the target and the sensor [m]
 * @param rho -> effective reflectance
 * @param area -> effective reflective area [m^2]
 * @return
 */
double computeIrradianceOfSatellite(const double r, const double rho, const double area);

double irradianceToApparentMagnitude(const double e);

class DistanceUpdateCallback : public osg::Drawable::UpdateCallback
{
public:
    DistanceUpdateCallback(osg::Camera* camera, osg::Light* light, osg::Geometry* geo, osg::Texture2D* tex);
    ~DistanceUpdateCallback();
    virtual void update(osg::NodeVisitor* nv, osg::Drawable*);
private:
    float computeCosineOfAngle(const osg::Vec3f& a, const osg::Vec3f& b);
    osg::Vec3d getCameraWorldPosition(osg::Camera* camera);
    void printAngle(const float cosin);
    double computeERA();
    void updateDarkPoint(double distance);
protected:
    osg::ref_ptr<osg::Camera> m_camera;
    osg::ref_ptr<osg::Light> m_light;
    osg::ref_ptr<osg::Geometry> m_geo;
    osg::ref_ptr<osg::Texture2D> m_tex;
    osg::ref_ptr<osg::Image> m_psfImage;    // the origional PSF image
    osg::Matrixd m_lastViewMatrix;
};

class TogglePBREventHandler : public osgGA::GUIEventHandler {
public:
    TogglePBREventHandler(osg::StateSet* stateset);

    virtual bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter&);

private:
    void applyPBR();

    void applyNonPBR();

    osg::ref_ptr<osg::StateSet> _stateset;
    bool _isPBR;
};

//osg::Vec4 computePixelSizeVector(const osg::Viewport& W, const osg::Matrix& P, const osg::Matrix& M)
//{
//    // pre adjust P00,P20,P23,P33 by multiplying them by the viewport window matrix.
//    // here we do it in short hand with the knowledge of how the window matrix is formed
//    // note P23,P33 are multiplied by an implicit 1 which would come from the window matrix.
//    // Robert Osfield, June 2002.
//
//    // scaling for horizontal pixels
//    float P00 = P(0, 0) * W.width() * 0.5f;
//    float P20_00 = P(2, 0) * W.width() * 0.5f + P(2, 3) * W.width() * 0.5f;
//    osg::Vec3 scale_00(M(0, 0) * P00 + M(0, 2) * P20_00,
//        M(1, 0) * P00 + M(1, 2) * P20_00,
//        M(2, 0) * P00 + M(2, 2) * P20_00);
//
//    // scaling for vertical pixels
//    float P10 = P(1, 1) * W.height() * 0.5f;
//    float P20_10 = P(2, 1) * W.height() * 0.5f + P(2, 3) * W.height() * 0.5f;
//    osg::Vec3 scale_10(M(0, 1) * P10 + M(0, 2) * P20_10,
//        M(1, 1) * P10 + M(1, 2) * P20_10,
//        M(2, 1) * P10 + M(2, 2) * P20_10);
//
//    float P23 = P(2, 3);
//    float P33 = P(3, 3);
//    osg::Vec4 pixelSizeVector(M(0, 2) * P23,
//        M(1, 2) * P23,
//        M(2, 2) * P23,
//        M(3, 2) * P23 + M(3, 3) * P33);
//
//    float scaleRatio = 0.7071067811f / sqrtf(scale_00.length2() + scale_10.length2());
//    pixelSizeVector *= scaleRatio;
//
//    return pixelSizeVector;
//}
}