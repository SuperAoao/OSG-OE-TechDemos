#version 120
varying vec2 TexCoords;
varying vec3 VetexPos;
varying vec3 Normal;
varying vec3 ViewDir;
varying mat3 TBNMatrix;
varying mat3 NormalWorldToModelMatrix;
varying vec3 ModelVetexPos;
varying vec3 ModelViewDir;
uniform samplerCube irradianceMap;
uniform samplerCube prefilterMap0;
uniform samplerCube prefilterMap1;
uniform samplerCube prefilterMap2;
uniform samplerCube prefilterMap3;
uniform samplerCube prefilterMap4;
uniform sampler2D brdfLUT;
uniform sampler2D albedoMap;
uniform sampler2D normalMap;
uniform sampler2D metallicMap;
uniform sampler2D roughnessMap;
uniform sampler2D aoMap;
uniform bool albedoSwitch;
uniform bool normalSwitch;
uniform bool metallicSwitch;
uniform bool roughnessSwitch;
uniform bool aoSwitch;
uniform bool albedoSRGBSwitch;
uniform bool IBLSwitch;
uniform bool gammaSwitch;
uniform bool hdrSwitch;
uniform bool pbrSwitch;
uniform bool localSkySwitch;
uniform vec3 albedoValue;
uniform vec3 normalValue;
uniform float metallicValue;
uniform float roughnessValue;
uniform float aoValue;
#define MAXLIGHTNUM 200
uniform vec4 globalLightPositions[MAXLIGHTNUM]; 
uniform vec3 globalLightColors[MAXLIGHTNUM];
uniform int globalLightNum;
uniform vec4 localLightPositions[MAXLIGHTNUM]; 
uniform vec3 localLightColors[MAXLIGHTNUM];
uniform int localLightNum;
uniform mat4 osg_ViewMatrix;
uniform mat4 osg_ViewMatrixInverse;
uniform vec3 envLighting;
const float PI = 3.14159265359;
vec3 worldNormal;
vec3 modelNormal;
mat3 mat3Func(mat4 m)
{
	mat3 mat3Matrix;
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			mat3Matrix[i][j] = m[i][j];
		}
	}
	return mat3Matrix;
}
vec3 getNormalFromMap()
{
	vec3 tangentNormal = normalValue;
	if(normalSwitch && pbrSwitch)
	{
		tangentNormal = texture2D(normalMap, TexCoords).xyz * 2.0 - 1.0;
	}
	worldNormal = normalize(TBNMatrix * tangentNormal);
	modelNormal = NormalWorldToModelMatrix * worldNormal;
	return mat3Func(osg_ViewMatrix) * worldNormal;
}
float DistributionGGX(vec3 N, vec3 H, float roughness)
{
	float a = roughness*roughness;
	float a2 = a*a;
	float NdotH = max(dot(N, H), 0.0);
	float NdotH2 = NdotH*NdotH;
	float nom   = a2;
	float denom = (NdotH2 * (a2 - 1.0) + 1.0);
	denom = PI * denom * denom;
	return nom / denom;
}
float GeometrySchlickGGX(float NdotV, float roughness)
{
	float r = (roughness + 1.0);
	float k = (r*r) / 8.0;
	float nom   = NdotV;
	float denom = NdotV * (1.0 - k) + k;
	return nom / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
	float NdotV = max(dot(N, V), 0.0);
	float NdotL = max(dot(N, L), 0.0);
	float ggx2 = GeometrySchlickGGX(NdotV, roughness);
	float ggx1 = GeometrySchlickGGX(NdotL, roughness);
	return ggx1 * ggx2;
}
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
	return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}
vec3 fresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness)
{
	return F0 + (max(vec3(1.0 - roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
}
vec3 prefilterMapLod(vec3 uvw, float roughnessLod)
{
	//return textureCube(prefilterMap2, uvw).xyz;
	vec3 startColor;
	vec3 endColor;
	float lerpW;
	if(roughnessLod <= 0.0)
	{
		return textureCube(prefilterMap0, uvw).xyz;
	}
	else if(roughnessLod <= 1.0)
	{
		startColor = textureCube(prefilterMap0, uvw).xyz + envLighting;
		endColor = textureCube(prefilterMap1, uvw).xyz + envLighting;
		lerpW = roughnessLod;
	}
	else if(roughnessLod <= 2.0)
	{
		startColor = textureCube(prefilterMap1, uvw).xyz + envLighting;
		endColor = textureCube(prefilterMap2, uvw).xyz + envLighting;
		lerpW = roughnessLod-1.0;
	}
	else if(roughnessLod <= 3.0)
	{
		startColor = textureCube(prefilterMap2, uvw).xyz + envLighting;
		endColor = textureCube(prefilterMap3, uvw).xyz + envLighting;
		lerpW = roughnessLod-2.0;
	}
	else
	{
		startColor = textureCube(prefilterMap3, uvw).xyz + envLighting;
		endColor = textureCube(prefilterMap4, uvw).xyz + envLighting;
		lerpW = roughnessLod-3.0;
	}

	return startColor*(1.0-lerpW)+endColor*lerpW;
}
void main()
{
	vec3 albedo = albedoValue;
	if(albedoSwitch)
	{
		albedo = texture2D(albedoMap, TexCoords).xyz;
		if(!albedoSRGBSwitch && pbrSwitch)
		{
			albedo = pow(albedo, vec3(2.2));
		}
	}
	vec3 N = getNormalFromMap();
	vec3 V = normalize(ViewDir);
	vec3 R = reflect(-V, N);
	if(!pbrSwitch)
	{
		vec3 ambient = vec3(0.2);
		vec3 diffuse = vec3(0.0);
		vec3 specular = vec3(0.0);
		for(int i = 0; i < globalLightNum; ++i)
		{
			vec4 lightViewPos = osg_ViewMatrix * globalLightPositions[i];
			vec3 L = normalize(lightViewPos.xyz);
			lightViewPos.w = globalLightPositions[i].w;
			float attenuation = 1.0;
			if(lightViewPos.w > 0.0)
			{
				L = normalize(lightViewPos.xyz - VetexPos);
				float distance = length(globalLightPositions[i].xyz - VetexPos);
				attenuation = 1.0 / (distance * distance);
			}
			float diff = max(dot(N, L), 0.0);
			vec3 reflectDir = reflect(-L, N);
			float spec = pow(max(dot(V, reflectDir), 0.0), 10.0);
			diffuse = diffuse + vec3(diff) * vec3(0.4);
			specular = specular + vec3(spec) * vec3(0.5);
		}
		gl_FragColor = vec4(albedo * (ambient + diffuse + specular), 1.0);
		return;
	}
	float metallic = metallicValue;
	if(metallicSwitch)
	{
		metallic = texture2D(metallicMap, TexCoords).x;
	}
	float roughness = roughnessValue;
	if(roughnessSwitch)	{
		roughness = texture2D(roughnessMap, TexCoords).x;
	}
	float ao = aoValue;
	if(aoSwitch)	{
		ao = texture2D(aoMap, TexCoords).x;
	}
	vec3 F0 = vec3(0.04);
	F0 = mix(F0, albedo, metallic);
	vec3 Lo = vec3(0.0);
	for(int i = 0; i < globalLightNum; ++i)
	{
		vec4 lightViewPos = osg_ViewMatrix * globalLightPositions[i];
		vec3 L = normalize(lightViewPos.xyz);
		lightViewPos.w = globalLightPositions[i].w;
		float attenuation = 1.0;
		if(lightViewPos.w > 0.0)
		{
			L = normalize(lightViewPos.xyz - VetexPos);
		}
		vec3 H = normalize(V + L);
		vec3 radiance = globalLightColors[i] * attenuation;
		float NDF = DistributionGGX(N, H, roughness);
		float G   = GeometrySmith(N, V, L, roughness);
		vec3 F    = fresnelSchlick(max(dot(H, V), 0.0), F0);
		vec3 nominator    = NDF * G * F;
		float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.001;
		vec3 specular = nominator / denominator;
		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;
		kD *= 1.0 - metallic;
		float NdotL = max(dot(N, L), 0.0);
		Lo += (kD * albedo / PI + specular) * radiance * NdotL;
	}
	vec3 modelN = normalize(modelNormal);
	vec3 modelV = normalize(ModelViewDir);
	vec3 modelR = reflect(-modelV, modelN);
	for(int i = 0; i < localLightNum; ++i)
	{
		vec4 lightViewPos = localLightPositions[i];
		vec3 modelL = normalize(lightViewPos.xyz);
		lightViewPos.w = localLightPositions[i].w;
		float attenuation = 1.0;
		if(lightViewPos.w > 0.0)
		{
			modelL = normalize(lightViewPos.xyz - ModelVetexPos);
			float distance = length(localLightPositions[i].xyz - ModelVetexPos);
			attenuation = 1.0 / (distance * distance);
		}
		vec3 modelH = normalize(modelV + modelL);
		vec3 radiance = localLightColors[i] * attenuation;
		float NDF = DistributionGGX(modelN, modelH, roughness);
		float G   = GeometrySmith(modelN, modelV, modelL, roughness);
		vec3 F    = fresnelSchlick(max(dot(modelH, modelV), 0.0), F0);
		vec3 nominator    = NDF * G * F;
		float denominator = 4.0 * max(dot(modelN, modelV), 0.0) * max(dot(modelN, modelL), 0.0) + 0.001;
		vec3 specular = nominator / denominator;
		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;
		kD *= 1.0 - metallic;
		float NdotL = max(dot(modelN, modelL), 0.0);
		Lo += (kD * albedo / PI + specular) * radiance * NdotL;
	}
	vec3 F = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, roughness);
	vec3 kS = F;
	vec3 kD = 1.0 - kS;
	kD *= 1.0 - metallic;
	vec3 ambient = vec3(0.0, 0.0, 0.0);
	if(IBLSwitch)	{
		vec3 IBLNormal;
		vec3 IBLR;
		if(localSkySwitch)
		{
			IBLNormal = modelNormal;
			IBLR = modelR;
		}
		else
		{
			IBLNormal = worldNormal ;
			IBLR = mat3Func(osg_ViewMatrixInverse) * R;
		}
		vec3 irradiance = textureCube(irradianceMap, IBLNormal).xyz + envLighting ;
		vec3 diffuse      = irradiance * albedo;
		const float MAX_REFLECTION_LOD = 4.0;
		vec3 prefilteredColor = prefilterMapLod(IBLR,  roughness * MAX_REFLECTION_LOD);
		vec2 brdf  = texture2D(brdfLUT, vec2(max(dot(N, V), 0.0), roughness)).xy;
		vec3 specular = prefilteredColor * (F * brdf.x + brdf.y);
		ambient = (kD * diffuse + specular) * ao;
	}
	vec3 color = ambient + Lo;
	if(hdrSwitch){
		color = color / (color + vec3(1.0));
	};	if(gammaSwitch){
		color = pow(color, vec3(1.0 / 2.2));
	};	gl_FragColor = vec4(color , 1.0);
}
