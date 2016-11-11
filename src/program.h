#pragma once

#include <GL/glew.h>
#include <glm/mat4x4.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

class Program
{
public:
	Program(std::string vertexShaderFilename, std::string fragmentShaderFilename);
	~Program();

	void useProgram();
	void stopUseProgram();

	void loadModelViewMatrix(const glm::mat4 &modelView);
	void loadProjMatrix(glm::mat4 &proj);
	void loadColorUniform(glm::vec3 &color);

	void destroy();

private:
	GLuint m_program;
	void loadShaders(std::string vertexShaderFilename, std::string fragmentShaderFilename);
};

