#pragma once

#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osg/Texture2D>
#include <osg/TextureCubeMap>

/**
  * @class PBRStateSet
  * @brief 
  * @note 
  * @author 
*/
class  PBRStateSet : public osg::StateSet
{
public:
	PBRStateSet();
	~PBRStateSet();
	enum texFormat
	{
		PNG = 0,
		JPG = 1
	};
	// AlbedoMapFileName: "Base_color.*"
	// NormalMapFileName: "Normal.*"
	// MetallicMapFileName: "Metallic.*"
	// RoughnessMapFileName: "Roughness.*"
	// AOMapFileName: "AO.*"
	void loadTextures(const std::string& directoryPath, const texFormat texFormat = PNG);
	bool setAlbedoMap(const std::string& filePath, bool sRGB = false, bool bUnrefAfterApply = true);
	bool setNormalMap(const std::string& filePath, bool sRGB = false, bool bUnrefAfterApply = true);
	bool setMetallicMap(const std::string& filePath, bool sRGB = false, bool bUnrefAfterApply = true);
	bool setRoughnessMap(const std::string& filePath, bool sRGB = false, bool bUnrefAfterApply = true);
	bool setAOMap(const std::string& filePath, bool sRGB = false, bool bUnrefAfterApply = true);
	static bool setEnvironmentLighting(osg::Vec3f lightPos);
private:
	// init shader program
	bool initProgram();
	osg::Texture2D* createTexture2D(osg::Image* image);

	// PBR tex
	osg::ref_ptr<osg::Texture2D> m_rpAlbedoMap;
	osg::ref_ptr<osg::Texture2D> m_rpAoMap;
	osg::ref_ptr<osg::Texture2D> m_rpMetallicMap;
	osg::ref_ptr<osg::Texture2D> m_rpNormalMap;
	osg::ref_ptr<osg::Texture2D> m_rpRoughnessMap;

	//PBR tex default val
	osg::ref_ptr<osg::Uniform> m_rpAlbedoValue;
	osg::ref_ptr<osg::Uniform> m_rpAoValue;
	osg::ref_ptr<osg::Uniform> m_rpMetallicValue;
	osg::ref_ptr<osg::Uniform> m_rpNormalValue;
	osg::ref_ptr<osg::Uniform> m_rpRoughnessValue;

	// global uniform
	static osg::ref_ptr<osg::Uniform> s_rpLightNum;
	static osg::ref_ptr<osg::Uniform> s_rpLightPos;
	static osg::ref_ptr<osg::Uniform> s_rpLightColors;
	static osg::ref_ptr<osg::Uniform> s_rpEnvLighting; 

	// switcher
	osg::ref_ptr<osg::Uniform> m_rpAlbedoSwitch;
	osg::ref_ptr<osg::Uniform> m_rpAoSwitch;
	osg::ref_ptr<osg::Uniform> m_rpMetallicSwitch;
	osg::ref_ptr<osg::Uniform> m_rpNormalSwitch;
	osg::ref_ptr<osg::Uniform> m_rpRoughnessSwitch;
	osg::ref_ptr<osg::Uniform> m_rpAlbedoSRGBSwitch;
	osg::ref_ptr<osg::Uniform> m_rpAutoRotateSkySwitch;

	// global switcher
	static osg::ref_ptr<osg::Uniform> s_rpGammaSwitch;
	static osg::ref_ptr<osg::Uniform> s_rpHdrSwitch;
	static osg::ref_ptr<osg::Uniform> s_rpPbrSwitch;
	static osg::ref_ptr<osg::Uniform> s_rpIBLSwitch;

	// IBL(need preCal)
	static osg::ref_ptr<osg::TextureCubeMap> s_rpEnvCubemap;
	static osg::ref_ptr<osg::TextureCubeMap> s_rpIrrCubemap;
	static osg::ref_ptr<osg::TextureCubeMap> s_rpPrefilterCubemap[5];
	static osg::ref_ptr<osg::Texture2D> s_rpBrdfLutTexture;

	static osg::Program* s_rpProgram;
	static bool	s_bIniProgram;	// is shader program initiliazed
	static bool	s_bCoreProfile;	// core			 -> glsl version 330
								// compatibility -> glsl version 120
};
