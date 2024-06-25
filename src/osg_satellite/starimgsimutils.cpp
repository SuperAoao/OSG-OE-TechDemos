#include "starimgsimutils.h"

#include <osgDB/ReadFile>
#include <osgDB/WriteFile>

namespace StarImgSimUtils
{
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

    DistanceUpdateCallback::DistanceUpdateCallback(osg::Camera* camera, osg::Light* light, osg::Geometry* geo, osg::Texture2D* tex)
        :m_camera(camera), m_light(light), m_geo(geo), m_tex(tex), m_psfImage(nullptr)
    {
        osg::Object* imageObj = m_tex->getImage()->clone(osg::CopyOp::DEEP_COPY_IMAGES);
        osg::Image* image = dynamic_cast<osg::Image*>(imageObj);
        if (image)
        {
            m_psfImage = image;
        }
        m_lastViewMatrix = m_camera->getViewMatrix();
    }
    DistanceUpdateCallback::~DistanceUpdateCallback()
    {
    }
    void DistanceUpdateCallback::update(osg::NodeVisitor* nv, osg::Drawable*)
    {
        osg::Matrixd curMatrix = m_camera->getViewMatrix();
        if (m_lastViewMatrix != curMatrix)
        {
            std::cout << "Moving!" << std::endl;
            double distance = curMatrix.getTrans().length();
            std::cout << "distance: " << distance << std::endl;
            if (distance < 10000)
            {
                double max = 7438980;
                double min = 3505390;
                double range = max - min;
                double factor = range / 7000;
                distance = min + factor * (10000 - distance);
            }
            updateDarkPoint(distance);
        }
        m_lastViewMatrix = curMatrix;

    }
    float DistanceUpdateCallback::computeCosineOfAngle(const osg::Vec3f& a, const osg::Vec3f& b)
    {
        osg::Vec3f vectorA = a;
        osg::Vec3f vectorB = b;

        // Normalize the vectors
        vectorA.normalize();
        vectorB.normalize();

        // Calculate the dot product
        float dotProduct = vectorA * vectorB; // Using operator* for dot product

        // Clamp the dot product within the valid range for acos

        dotProduct = osg::clampTo(dotProduct, -1.0f, 1.0f);
        return dotProduct;
    }
    osg::Vec3d DistanceUpdateCallback::getCameraWorldPosition(osg::Camera* camera)
    {
        osg::Matrixd viewMatrix = camera->getViewMatrix();
        osg::Matrixd inverseViewMatrix;
        inverseViewMatrix.invert(viewMatrix);

        // Extract translation component from the inverted view matrix
        osg::Vec3d cameraWorldPosition = inverseViewMatrix.getTrans();

        return cameraWorldPosition;
    }
    void DistanceUpdateCallback::printAngle(const float cosin)
    {
        // Calculate the angle in radians
        float angleRadians = std::acos(cosin);

        // Convert the angle to degrees (optional)
        float angleDegrees = osg::RadiansToDegrees(angleRadians);

        // Output the results
        std::cout << "Angle between vectors: " << angleRadians << " radians (" << angleDegrees << " degrees)" << std::endl;
    }
    double DistanceUpdateCallback::computeERA()
    {
        double era = 0.0;
        // pos of sun
        osg::Vec3d sunPos = osg::Vec3d(m_light->getPosition().x(), m_light->getPosition().y(), m_light->getPosition().z());
        // pos of cam
        osg::Vec3d camPos = getCameraWorldPosition(m_camera);
        std::cout << "current Cam Pos: (" << camPos.x() << ", " << camPos.y() << ", " << camPos.z() << ")" << std::endl;
        camPos.normalize();
        std::vector<osg::Vec3d> vec;
        osg::Vec3d topN(0.0, 0.0, 1.0);
        osg::Vec3d bottomN(0.0, 0.0, -1.0);
        osg::Vec3d leftN(-1.0, 0.0, 0.0);
        osg::Vec3d rightN(1.0, 0.0, 0.0);
        osg::Vec3d frontN(0.0, -1.0, 0.0);
        osg::Vec3d backN(0.0, 1.0, 0.0);
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

            std::cout << "current Sun Pos: " << sunPos.x() << ", " << sunPos.y() << ", " << sunPos.z() << std::endl;
            std::cout << "currentNormal:( " << normal.x() << ", " << normal.y() << ", " << normal.z() << ")" << std::endl;
            double cos1i = computeCosineOfAngle(camPos, normal);
            std::cout << "cos1i: " << cos1i << std::endl;
            printAngle(cos1i);
            double cos2i = computeCosineOfAngle(osg::Vec3d(sunPos.x(), sunPos.y(), sunPos.z()), normal);
            std::cout << "cos2i: " << cos2i << std::endl;
            printAngle(cos2i);
            double area = planeArea * cos1i * cos2i;
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
    void DistanceUpdateCallback::updateDarkPoint(double distance)
    {
        double E = computeIrradianceOfSatellite(distance, g_rho, computeERA());
        double target_mag = irradianceToApparentMagnitude(E);
        double d_len = 0.03; // unit m
        double t_exposure = 5;    // unit s
        double n_pe = getPhotonsFromMagnitude(target_mag, d_len, t_exposure);
        // apply energy
        double grayRatio = n_pe / StarImgSimUtils::g_n_e;
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
    TogglePBREventHandler::TogglePBREventHandler(osg::StateSet* stateset)
        : _stateset(stateset), _isPBR(true)
    {
        // Initial state is PBR
        applyPBR();
    }
    bool TogglePBREventHandler::handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter&)
    {
        if (ea.getEventType() == osgGA::GUIEventAdapter::KEYDOWN) {
            if (ea.getKey() == 'f' || ea.getKey() == 'F') {
                _isPBR = !_isPBR;
                if (_isPBR) {
                    applyPBR();
                }
                else {
                    applyNonPBR();
                }
                return true; // event handled
            }
        }
        return false; // event not handled
    }
    void TogglePBREventHandler::applyPBR()
    {
    }
    void TogglePBREventHandler::applyNonPBR()
    {
    }
}