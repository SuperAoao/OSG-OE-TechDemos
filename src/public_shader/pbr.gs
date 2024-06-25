#version 120
#extension GL_EXT_geometry_shader4 : enable
varying in vec2 vTexCoords[];
varying out vec2 TexCoords;
varying in vec3 vVetexPos[];
varying out vec3 VetexPos;
varying in vec3 vNormal[];
varying out vec3 Normal;
varying in vec3 vViewDir[];
varying out vec3 ViewDir;
varying out mat3 TBNMatrix;
varying out mat3 NormalWorldToModelMatrix;
varying in mat3 vNormalModelToWorldMatrix[];
varying in mat3 vNormalWorldToModelMatrix[];
varying in vec3 vModelViewDir[];
varying out vec3 ModelViewDir;
varying in vec3 vModelVetexPos[];
varying out vec3 ModelVetexPos;
void main(void)
{
	vec3 edge1 = gl_PositionIn[1].xyz - gl_PositionIn[0].xyz;
	vec3 edge2 = gl_PositionIn[2].xyz - gl_PositionIn[0].xyz;
	vec2 deltaUV1 = vTexCoords[1] - vTexCoords[0];
	vec2 deltaUV2 = vTexCoords[2] - vTexCoords[0];
	float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
;	vec3 tangent;
	tangent.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
;	tangent.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
;	tangent.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);
;	tangent = normalize(tangent);
;    for(int i=0;i< 3; i++)
    {
        TexCoords = vTexCoords[i];
        VetexPos = vVetexPos[i];
        Normal = vNormal[i];
        ViewDir = vViewDir[i];
        ModelVetexPos = vModelVetexPos[i];
        ModelViewDir = vModelViewDir[i];
        vec3 B  = -normalize(cross(Normal, tangent));
        TBNMatrix = mat3(vNormalModelToWorldMatrix[i] * tangent, vNormalModelToWorldMatrix[i] * B, vNormalModelToWorldMatrix[i] * Normal);
        NormalWorldToModelMatrix = vNormalWorldToModelMatrix[i];
        gl_Position = gl_ModelViewProjectionMatrix * gl_PositionIn[i];
        EmitVertex();    }
    EndPrimitive();
}
