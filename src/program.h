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

	void useProgram() const;
	void stopUseProgram() const;

	void loadModelViewMatrix(const glm::mat4 &modelView) const;
	void loadProjMatrix(const glm::mat4 &proj) const;
	void loadColorUniform(const glm::vec4 &color) const;

	void destroy();

private:
	GLuint m_program;
	void loadShaders(std::string vertexShaderFilename, std::string fragmentShaderFilename);
};

