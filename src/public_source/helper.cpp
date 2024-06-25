#include "../public_header/helper.h"

 std::string Helper::readShaderFile(const std::string& filePath) {
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