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
#include <osgDB/WriteFile>
#include <osg/Image>
#include <osg/NodeVisitor>
#include <osg/Drawable>
#include <osg/LineWidth>
#include <osgViewer/ViewerEventHandlers>

// exposure time 1s, d_length 1 mm^2
const double g_CCDphotons = 19100;
// satellite
const double g_E_sun = 573; // Irradiance of the sun measured in[W/m^2]
const double g_sun_app_mag = -26.7; // sun apparent magnitude
const double g_rho = 0.25; // effective reflectance
// CCD gray parameters
const double n_ec = 45000;    // CCD per pixel max photons
const double n_e = n_ec / 256;    // photons per gray scale


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

// we consider center point (xc,yc) = (0,0)
double simplifiedGaussianFunc(const double x, const double y, const double xc, const double yc, const double sigma)
{
    double xDis = x - xc;
    double yDis = y - yc;
    return std::exp(-(xDis * xDis + yDis * yDis) / (2 * sigma * sigma));
}

double getEnergyProportion(const double x_start, const double y_start, const double x_end, const double y_end, const double steps, const double sigma)
{
    // gaussian spread center coords are (0,0)
    double hx = (x_end - x_start) / steps; // x 方向步长
    double hy = (y_end - y_start) / steps; // y 方向步长
    double xc = 0.0;
    double yc = 0.0;
    double sum = 0.0;
    for (uint32_t i = 0; i < steps; ++i)
    {
        for (uint32_t j = 0; j < steps; ++j)
        {
            double x0 = x_start + i * hx; // 当前子区间左边界的 x 值
            double x1 = x_start + (i + 1) * hx; // 当前子区间右边界的 x 值
            double y0 = y_start + j * hy; // 当前子区间下边界的 y 值
            double y1 = y_start + (j + 1) * hy; // 当前子区间上边界的 y 值

            sum += (simplifiedGaussianFunc(x0, y0, xc, yc, sigma) +
                simplifiedGaussianFunc(x1, y0, xc, yc, sigma) +
                simplifiedGaussianFunc(x0, y1, xc, yc, sigma) +
                simplifiedGaussianFunc(x1, y1, xc, yc, sigma)) / 4.0 * hx * hy;
        }
    }
    sum = 1 / (2 * osg::PI * sigma * sigma) * sum;
    return sum;
}

osg::Texture2D* genPSFTexture(unsigned int pixelSize, const double sigma)
{
    osg::Image* image = new osg::Image;
    image->allocateImage(pixelSize, pixelSize, 1, GL_LUMINANCE, GL_UNSIGNED_BYTE);
    unsigned char* data = image->data();
    double sum = 0.0;
    double centerIndex = pixelSize / 2;
    for (unsigned int i = 0; i < pixelSize; i++)
    {
        for (unsigned int j = 0; j < pixelSize; j++)
        {
            double energyProportion = getEnergyProportion((double)i - centerIndex, (double)j - centerIndex, 
                (double)i - centerIndex + 1, (double)j - centerIndex + 1, 100, sigma);
            sum += energyProportion;
            //energyProportion = 182745.87074774309 * energyProportion / n_e / 256.0;
            //energyProportion = 487315350.22900116 * energyProportion / n_e / 256.0;
            data[i * pixelSize + j] = (unsigned char)(255.99f * std::min(energyProportion, 1.0));
        }
    }

    osg::Texture2D* texture = new osg::Texture2D;
    texture->setImage(image);
    texture->setInternalFormat(GL_RGB);
    texture->setBorderColor(osg::Vec4(0.0, 0.0, 0.0, 0.0));
    texture->setFilter(osg::Texture2D::MIN_FILTER, osg::Texture2D::LINEAR_MIPMAP_LINEAR);
    texture->setFilter(osg::Texture2D::MAG_FILTER, osg::Texture2D::LINEAR);
    texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::CLAMP_TO_BORDER);
    texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::CLAMP_TO_BORDER);

    osgDB::writeImageFile(*image, std::string("test.png"));
    return texture;
}



// get photons from star magitude on CCD
double getPhotonsFromMagnitude(double mv, double d_len, double t_in)
{
    double n_pe = g_CCDphotons * 1 / std::pow(2.5, mv) * t_in * osg::PI * (d_len / 2) * (d_len / 2) * 1e6;
    return n_pe;
}

/**
 * @brief
 * @param r -> the distance between the target and the sensor [m]
 * @param rho -> effective reflectance
 * @param area -> effective reflective area [m^2]
 * @return
 */
double computeIrradianceOfSatellite(const double r, const double rho, const double area)
{
    double E = g_E_sun / (osg::PI * r * r) * rho * area;
    return E;
}

double irradianceToApparentMagnitude(const double e)
{
    double m = 2.5 * std::log10(g_E_sun / e) + g_sun_app_mag;
    return m;
}

class DistanceUpdateCallback : public osg::Drawable::UpdateCallback
{
public:
    DistanceUpdateCallback(osg::Camera* camera, osg::Geometry* geo, osg::Texture2D* tex)
        : m_camera(camera), m_geo(geo) , m_tex(tex), m_psfImage(nullptr)
    {
        osg::Object* imageObj = m_tex->getImage()->clone(osg::CopyOp::DEEP_COPY_IMAGES);
        osg::Image* image = dynamic_cast<osg::Image*>(imageObj);
        if (image)
        {
            m_psfImage = image;
        }
        m_lastViewMatrix = m_camera->getViewMatrix();
    }

    virtual void update(osg::NodeVisitor* nv, osg::Drawable*)
    {
        osg::Matrixd curMatrix = m_camera->getViewMatrix();
        if (m_lastViewMatrix != curMatrix)
        {
            std::cout << "Moving!" << std::endl;
            double distance = curMatrix.getTrans().length();
            std::cout << "distance: " << distance << std::endl;
            updateDarkPoint(distance);
        }
        m_lastViewMatrix = curMatrix;
        
    }
private:
    double computeCosineOfAngle(const osg::Vec3& a, const osg::Vec3& b)
    {
        // Compute dot product of the vectors
        double dotProduct = a * b;

        // Compute magnitudes of the vectors
        double magnitudeA = a.length();
        double magnitudeB = b.length();

        // Compute cosine of the angle between the vectors
        double cosine = dotProduct / (magnitudeA * magnitudeB);

        return cosine;
    }
    osg::Vec3d getCameraWorldPosition(osg::Camera* camera)
    {
        osg::Matrixd viewMatrix = camera->getViewMatrix();
        osg::Matrixd inverseViewMatrix;
        inverseViewMatrix.invert(viewMatrix);

        // Extract translation component from the inverted view matrix
        osg::Vec3d cameraWorldPosition = inverseViewMatrix.getTrans();

        return cameraWorldPosition;
    }
    double computeERA()
    {
        double era = 0.0;
        // pos of sun
        osg::Vec3d sunPos(0.0, 1.0, 0.0);
        // pos of cam
        osg::Vec3d camPos = getCameraWorldPosition(m_camera);
        std::cout << "current Cam Pos: (" << camPos.x() << ", " << camPos.y() << ", " << camPos.z() << ")" << std::endl;
        camPos.normalize();
        std::vector<osg::Vec3d> vec;
        osg::Vec3d topN(0.0, 1.0, 0.0);
        osg::Vec3d bottomN(0.0, -1.0, 0.0);
        osg::Vec3d leftN(-1.0, 0.0, 0.0);
        osg::Vec3d rightN(1.0, 0.0, 0.0);
        osg::Vec3d frontN(0.0, 0.0, 1.0);
        osg::Vec3d backN(0.0, 0.0, -1.0);
        vec.push_back(topN);
        vec.push_back(bottomN);
        vec.push_back(leftN);
        vec.push_back(rightN);
        vec.push_back(frontN);
        vec.push_back(backN);

        double planeArea = 1.0; // m^2
        double swingsArea = 5.0 * 2; // m^2
        
        for (int i = 0; i < vec.size(); ++i)
        {
            osg::Vec3d normal = vec.at(i);
            std::cout << "currentNormal:( " << normal.x() << ", " << normal.y() << ", " << normal.z() << ")" << std::endl;
            double cos1i = computeCosineOfAngle(camPos, normal);
            double cos2i = computeCosineOfAngle(sunPos, normal);
            std::cout << "cos1i: " << cos1i << ", cos2i: " << cos2i << std::endl;
            double area = planeArea *  cos1i * cos2i;
            std::cout << "area(" << i << "): " << area << std::endl;
            if (area > 0)
            {
                era += area;
            }
        }
        double area = planeArea * computeCosineOfAngle(camPos, topN) * computeCosineOfAngle(sunPos, topN);
        std::cout << "Swings area: " << area << std::endl;
        if (area > 0)
        {
            era += area;
        }
        std::cout << "ERA: " << era << std::endl;
        return era;

    }

    void updateDarkPoint(double distance)
    {
        computeERA();
        double E = computeIrradianceOfSatellite(distance, g_rho, 0.800);
        double target_mag = irradianceToApparentMagnitude(E);
        double d_len = 0.03; // unit m
        double t_exposure = 5;    // unit s
        double n_pe = getPhotonsFromMagnitude(target_mag, d_len, t_exposure);
        // apply energy
        double grayRatio = n_pe / n_e;
        if (m_psfImage->valid())
        {
            unsigned char* imgData = m_psfImage->data();
            unsigned char* texImgData = m_tex->getImage()->data();
            if (nullptr == imgData || nullptr == texImgData)
            {
                return;
            }
            for (int i = 0; i < 12; ++i)
            {
                for (int j = 0; j < 12; ++j)
                {
                    double proportion = (double)imgData[i * 12 + j] / 256.0;
                    double grayVal = std::min(proportion * grayRatio, 255.0);
                    texImgData[i * 12 + j] = (unsigned char)(grayVal);
                }
            }
            m_tex->dirtyTextureObject();
        }
    };
protected:
    osg::ref_ptr<osg::Camera> m_camera;
    osg::ref_ptr<osg::Geometry> m_geo;
    osg::ref_ptr<osg::Texture2D> m_tex;
    osg::ref_ptr<osg::Image> m_psfImage;    // the origional PSF image
    osg::Matrixd m_lastViewMatrix;
};

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
    double E = computeIrradianceOfSatellite(20440.0, g_rho, 0.800);
    double target_mag = irradianceToApparentMagnitude(E);
    double d_len = 0.03; // unit m
    double t_exposure = 5;    // unit s
    double n_pe = getPhotonsFromMagnitude(target_mag, d_len, t_exposure);

    osg::ref_ptr<osg::Texture2D> PSFTex = genPSFTexture(12, 0.85);
    // apply energy
    double grayRatio = n_pe / n_e;
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
    osg::ref_ptr<osg::MatrixTransform> mtScaleNode = new osg::MatrixTransform;
    mtScaleNode->setMatrix(osg::Matrix::scale(10, 10, 10));
    mtScaleNode->addChild(axesNode);

    osg::ref_ptr<osg::MatrixTransform> satellitePosMt = new osg::MatrixTransform;
    satellitePosMt->addChild(osgDB::readNodeFile("GPS.ive"));
    SatelliteGroup->addChild(satellitePosMt.get());
    SatelliteGroup->addChild(mtScaleNode.get());

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
    
    program->addShader(new osg::Shader(osg::Shader::VERTEX, readShaderFile("D:\\DevGitHub\\OSG-OE-TechDemos\\src\\osgsatellite\\satellitepoint.vs").c_str())); // Load the vertex shader
    program->addShader(new osg::Shader(osg::Shader::FRAGMENT, readShaderFile("D:\\DevGitHub\\OSG-OE-TechDemos\\src\\osgsatellite\\satellitepoint.fs").c_str())); // Load the fragment shader

    // Attach the shader program to the geometry
    geom->getOrCreateStateSet()->setAttributeAndModes(program.get());
    // Get the width of the window
    //unsigned int width = viewer->getCamera()->getViewport()->width();
    geom->getOrCreateStateSet()->setTextureAttributeAndModes(0, PSFTex.get(), osg::StateAttribute::ON);
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_width_half", 1280)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("window_height_half", 720)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("tex", 0)); // for temporary
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("n_pe", n_pe));
    geom->getOrCreateStateSet()->addUniform(new osg::Uniform("n_e", n_e));
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
    geom->setUpdateCallback(new DistanceUpdateCallback(viewer->getCamera(), geom.get(), PSFTex.get()));
    satellitePointGeode->addDrawable(geom.get());
    satellitePointGeode->getOrCreateStateSet()->setRenderBinDetails(1, "RenderBin");
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Depth(), osg::StateAttribute::OFF);
    satellitePointGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::Light(), osg::StateAttribute::OFF);
    satellitePosMt->addChild(satellitePointGeode.get());
    rootGroup->addChild(SatelliteGroup.get());
    viewer->setSceneData(rootGroup.get());
    viewer->addEventHandler(new osgViewer::StatsHandler());
    double fovy, aspectRatio, zNear, zFar;
    viewer->getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, zFar);
    viewer->getCamera()->setCullingMode(osg::CullSettings::NO_CULLING);
    viewer->getCamera()->setProjectionMatrixAsPerspective(fovy, aspectRatio, zNear, 100000);
    return viewer->run();
}

