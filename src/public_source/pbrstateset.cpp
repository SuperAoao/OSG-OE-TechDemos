#include "pbrstateset.h"
#include "helper.h"

#include <cstdlib>

#include <osg/NodeVisitor>
#include <osgUtil/CullVisitor>
#include <osgDB/WriteFile>
#include <osg/TexGen>
#include <osg/TexEnv>
#include <osg/TexMat>
#include <osg/ShapeDrawable>
#include <osg/Depth>
#include <osg/Notify>

#define MAXLIGHTNUM 200

const std::string g_vsFileName = "../shaders/pbr.vs";
const std::string g_gsFileName = "../shaders/pbr.gs";
const std::string g_fsFileName = "../shaders/pbr.fs";
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpGammaSwitch = new osg::Uniform("gammaSwitch", true);
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpHdrSwitch = new osg::Uniform("hdrSwitch", true);
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpPbrSwitch = new osg::Uniform("pbrSwitch", true);
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpIBLSwitch = new osg::Uniform("IBLSwitch", bool(false));
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpLightNum = new osg::Uniform("globalLightNum", int(0));
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpLightPos = new osg::Uniform(osg::Uniform::FLOAT_VEC4, "globalLightPositions", MAXLIGHTNUM);
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpLightColors = new osg::Uniform(osg::Uniform::FLOAT_VEC3, "globalLightColors", MAXLIGHTNUM);
osg::ref_ptr<osg::Uniform> PBRStateSet::s_rpEnvLighting = new osg::Uniform("envLighting", osg::Vec3f(0.0f, 0.0f, 0.0f));

osg::ref_ptr<osg::TextureCubeMap> PBRStateSet::s_rpEnvCubemap = nullptr;
osg::ref_ptr<osg::TextureCubeMap> PBRStateSet::s_rpIrrCubemap = nullptr; 
osg::ref_ptr<osg::TextureCubeMap> PBRStateSet::s_rpPrefilterCubemap[5] = { nullptr, nullptr, nullptr, nullptr ,nullptr };
osg::ref_ptr<osg::Texture2D> PBRStateSet::s_rpBrdfLutTexture = nullptr; 

osg::Program* PBRStateSet::s_rpProgram = new osg::Program;
bool PBRStateSet::s_bIniProgram = false;
bool PBRStateSet::s_bCoreProfile = false;
PBRStateSet::PBRStateSet() :
	m_rpAlbedoMap(nullptr), 
	m_rpNormalMap(nullptr), 
	m_rpMetallicMap(nullptr), 
	m_rpRoughnessMap(nullptr), 
	m_rpAoMap(nullptr),
	m_rpAlbedoValue(nullptr), 
	m_rpNormalValue(nullptr), 
	m_rpMetallicValue(nullptr), 
	m_rpRoughnessValue(nullptr), 
	m_rpAoValue(nullptr),
	m_rpAlbedoSwitch(nullptr), 
	m_rpNormalSwitch(nullptr), 
	m_rpMetallicSwitch(nullptr), 
	m_rpRoughnessSwitch(nullptr), 
	m_rpAoSwitch(nullptr),
	m_rpAlbedoSRGBSwitch(nullptr), 
	m_rpAutoRotateSkySwitch(nullptr)
{
	setTextureAttributeAndModes(1, s_rpPrefilterCubemap[0], osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(2, s_rpPrefilterCubemap[1], osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(3, s_rpPrefilterCubemap[2], osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(4, s_rpPrefilterCubemap[3], osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(5, s_rpPrefilterCubemap[4], osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(6, s_rpBrdfLutTexture, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(7, m_rpAlbedoMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(8, m_rpNormalMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(9, m_rpMetallicMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(10, m_rpRoughnessMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(11, m_rpAoMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
	setTextureAttributeAndModes(12, s_rpIrrCubemap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);

	m_rpAlbedoSwitch = new osg::Uniform("albedoSwitch", false);
	m_rpAoSwitch = new osg::Uniform("aoSwitch", false);
	m_rpMetallicSwitch = new osg::Uniform("metallicSwitch", false);
	m_rpNormalSwitch = new osg::Uniform("normalSwitch", false);
	m_rpRoughnessSwitch = new osg::Uniform("roughnessSwitch", false);

	m_rpAlbedoSRGBSwitch = new osg::Uniform("albedoSRGBSwitch", true);
	m_rpAutoRotateSkySwitch = new osg::Uniform("localSkySwitch", false);

	m_rpAlbedoValue = new osg::Uniform("albedoValue", osg::Vec3(1.0, 1.0, 1.0));
	m_rpAoValue = new osg::Uniform("aoValue", 1.0f);
	m_rpMetallicValue = new osg::Uniform("metallicValue", 0.0f);
	m_rpNormalValue = new osg::Uniform("normalValue", osg::Vec3(0.0, 0.0, 1.0));
	m_rpRoughnessValue = new osg::Uniform("roughnessValue", 0.0f);

	addUniform(new osg::Uniform("irradianceMap", int(12)));
	addUniform(new osg::Uniform("prefilterMap0", int(1)));
	addUniform(new osg::Uniform("prefilterMap1", int(2)));
	addUniform(new osg::Uniform("prefilterMap2", int(3)));
	addUniform(new osg::Uniform("prefilterMap3", int(4)));
	addUniform(new osg::Uniform("prefilterMap4", int(5)));
	addUniform(new osg::Uniform("brdfLUT", int(6)));
	addUniform(new osg::Uniform("albedoMap", int(7)));
	addUniform(new osg::Uniform("normalMap", int(8)));
	addUniform(new osg::Uniform("metallicMap", int(9)));
	addUniform(new osg::Uniform("roughnessMap", int(10)));
	addUniform(new osg::Uniform("aoMap", int(11)));

	addUniform(m_rpAlbedoSwitch);
	addUniform(m_rpAoSwitch);
	addUniform(m_rpMetallicSwitch);
	addUniform(m_rpNormalSwitch);
	addUniform(m_rpRoughnessSwitch);
	addUniform(m_rpAlbedoSRGBSwitch);
	addUniform(m_rpAutoRotateSkySwitch);

	addUniform(s_rpGammaSwitch);
	addUniform(s_rpPbrSwitch);
	addUniform(s_rpHdrSwitch);
	addUniform(s_rpIBLSwitch);

	addUniform(m_rpAlbedoValue);
	addUniform(m_rpAoValue);
	addUniform(m_rpMetallicValue);
	addUniform(m_rpNormalValue);
	addUniform(m_rpRoughnessValue);

	addUniform(s_rpLightNum);
	addUniform(s_rpLightPos);
	addUniform(s_rpLightColors);
	addUniform(s_rpEnvLighting);

	initProgram();

}

PBRStateSet::~PBRStateSet()
{
}

void PBRStateSet::loadTextures(const std::string& directoryPath, const texFormat texFormat)
{
	switch (texFormat)
	{
	case PNG:
	{
		setAlbedoMap(std::string(directoryPath).append("\\Base_color.png"));
		setMetallicMap(std::string(directoryPath).append("\\Metallic.png"));
		setNormalMap(std::string(directoryPath).append("\\Normal.png"));
		setRoughnessMap(std::string(directoryPath).append("\\Roughness.png"));
		setAOMap(std::string(directoryPath).append("\\AO.png"));
	}break;
	case JPG:
	{
		setAlbedoMap(std::string(directoryPath).append("\\Base_color.jpg"));
		setMetallicMap(std::string(directoryPath).append("\\Metallic.jpg"));
		setNormalMap(std::string(directoryPath).append("\\Normal.jpg"));
		setRoughnessMap(std::string(directoryPath).append("\\Roughness.jpg"));
		setAOMap(std::string(directoryPath).append("\\AO.jpg"));
	}break;
	default:
		break;
	}
	
}

bool PBRStateSet::setAlbedoMap(const std::string& filePath, bool sRGB , bool bUnrefAfterApply )
{
	osg::Image* albedoImage = osgDB::readImageFile(filePath);
	if (albedoImage)
	{
		osg::ref_ptr<osg::Texture2D> texture = createTexture2D(albedoImage);
		if (texture)
		{
			m_rpAlbedoSwitch->set(true);
			m_rpAlbedoSRGBSwitch->set(sRGB);
			m_rpAlbedoMap = texture;
			m_rpAlbedoMap->setUnRefImageDataAfterApply(bUnrefAfterApply);
			setTextureAttributeAndModes(7, m_rpAlbedoMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			return true;
		}
	}
	osg::notify(osg::WARN) << "PBRStateSet::setAlbedoMap failed!" << std::endl;
	return false;
}

bool PBRStateSet::setNormalMap(const std::string& filePath, bool sRGB, bool bUnrefAfterApply)
{
	osg::Image* normalImage = osgDB::readImageFile(filePath);
	if (normalImage)
	{
		osg::ref_ptr<osg::Texture2D> texture = createTexture2D(normalImage);
		if (texture)
		{
			m_rpNormalSwitch->set(true);
			m_rpNormalMap = texture;
			m_rpNormalMap->setUnRefImageDataAfterApply(bUnrefAfterApply);
			setTextureAttributeAndModes(8, m_rpNormalMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			return true;
		}
	}
	osg::notify(osg::WARN) << "PBRStateSet::setNormalMap failed!" << std::endl;
	return false;
}

bool PBRStateSet::setMetallicMap(const std::string& filePath, bool sRGB, bool bUnrefAfterApply)
{
	osg::Image* metallicImage = osgDB::readImageFile(filePath);
	if (metallicImage)
	{
		osg::ref_ptr<osg::Texture2D> texture = createTexture2D(metallicImage);
		if (texture)
		{
			m_rpMetallicSwitch->set(true);
			m_rpMetallicMap = texture;
			m_rpMetallicMap->setUnRefImageDataAfterApply(bUnrefAfterApply);
			setTextureAttributeAndModes(9, m_rpMetallicMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			return true;
		}
	}
	osg::notify(osg::WARN) << "PBRStateSet::setMetallicMap failed!" << std::endl;
	return false;
}

bool PBRStateSet::setRoughnessMap(const std::string& filePath, bool sRGB, bool bUnrefAfterApply)
{
	osg::Image* roughnessImage = osgDB::readImageFile(filePath);
	if (roughnessImage)
	{
		osg::ref_ptr<osg::Texture2D> texture = createTexture2D(roughnessImage);
		if (texture)
		{
			m_rpRoughnessSwitch->set(true);
			m_rpRoughnessMap = texture;
			m_rpRoughnessMap->setUnRefImageDataAfterApply(bUnrefAfterApply);
			setTextureAttributeAndModes(10, m_rpRoughnessMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			return true;
		}
	}
	return false;
}

bool PBRStateSet::setAOMap(const std::string& filePath, bool sRGB, bool bUnrefAfterApply)
{
	osg::Image* aoImage = osgDB::readImageFile(filePath);
	if (aoImage)
	{
		osg::ref_ptr<osg::Texture2D> texture = createTexture2D(aoImage);
		if (texture)
		{
			m_rpAoSwitch->set(true);
			m_rpAoMap = texture;
			m_rpAoMap->setUnRefImageDataAfterApply(bUnrefAfterApply);
			setTextureAttributeAndModes(11, m_rpAoMap, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			return true;
		}
	}
	return false;
}

bool PBRStateSet::setEnvironmentLighting(osg::Vec3f lightPos)
{
	if (s_rpEnvLighting)
	{
		s_rpEnvLighting->set(lightPos);
		return true;
	}
	return false;
}

bool PBRStateSet::initProgram()
{
	if (s_bIniProgram)
	{
		osg::notify(osg::DEBUG_INFO) << "PBR shader program has already been initialized!" << std::endl;
		return true;
	}
	else
	{
		if (s_bCoreProfile)
		{
			osg::notify(osg::FATAL) << "PBR shader program initialize failed: glsl330 version is not implemented" << std::endl;
			return false;
		}
		else
		{
			s_rpProgram->addShader(new osg::Shader(osg::Shader::VERTEX, Helper::readShaderFile(g_vsFileName)));
			s_rpProgram->addShader(new osg::Shader(osg::Shader::GEOMETRY, Helper::readShaderFile(g_gsFileName)));
			//设置几何着色器属性
			s_rpProgram->setParameter(GL_GEOMETRY_VERTICES_OUT_EXT, 3);
			s_rpProgram->setParameter(GL_GEOMETRY_INPUT_TYPE_EXT, GL_TRIANGLES);
			//program->setParameter( GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLES );

			s_rpProgram->addShader(new osg::Shader(osg::Shader::FRAGMENT, Helper::readShaderFile(g_fsFileName)));
			setAttributeAndModes(s_rpProgram, osg::StateAttribute::ON | osg::StateAttribute::PROTECTED);
			s_bIniProgram = true;
			return true;
		}
	}
}

osg::Texture2D* PBRStateSet::createTexture2D(osg::Image* image)
{
	osg::Texture2D* texture = nullptr;
	if (image)
	{
		texture = new osg::Texture2D;
		texture->setImage(image);
		//texture->setDataVariance(osg::Object::DYNAMIC);
		texture->setWrap(osg::Texture2D::WRAP_S, osg::Texture2D::REPEAT);
		texture->setWrap(osg::Texture2D::WRAP_T, osg::Texture2D::REPEAT);
		texture->setFilter(osg::Texture::MIN_FILTER, osg::Texture::LINEAR_MIPMAP_LINEAR);
		texture->setFilter(osg::Texture::MAG_FILTER, osg::Texture::LINEAR);
		//texture->setMaxAnisotropy(16.0f);
	}
	return texture;
}
