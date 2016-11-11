#include "program.h"

Program::Program(std::string vertexShaderFilename, std::string fragmentShaderFilename)
{
	m_program = glCreateProgram();
	loadShaders(vertexShaderFilename, fragmentShaderFilename);

}

Program::~Program()
{
}

void Program::useProgram()
{
	glUseProgram(m_program);
}

void Program::stopUseProgram()
{
	glUseProgram(0);
}

void Program::loadModelViewMatrix(const glm::mat4 & modelView)
{
	GLint matrixLocation = glGetUniformLocation(m_program, "modelView");
	if (matrixLocation < 0) {

		std::cout << "Invalid uniform" << std::endl;
	}
	glUniformMatrix4fv(matrixLocation, 1, GL_FALSE, glm::value_ptr(modelView));
}

void Program::loadProjMatrix(glm::mat4 & proj)
{
	GLint matrixLocation = glGetUniformLocation(m_program, "proj");
	if (matrixLocation < 0) {

		std::cout << "Invalid uniform" << std::endl;
	}
	glUniformMatrix4fv(matrixLocation, 1, GL_FALSE, glm::value_ptr(proj));
}

void Program::loadColorUniform(glm::vec3 & color)
{
	GLint colorLocation = glGetUniformLocation(m_program, "color");
	glUniform3f(colorLocation, color[0], color[1], color[2]);
}

void Program::destroy()
{
	glDeleteProgram(m_program);
}

void Program::loadShaders(std::string vertexShaderFilename, std::string fragmentShaderFilename)
{
	GLuint vertex_shader, fragment_shader;

	std::ifstream is(vertexShaderFilename);
	if (is) {
		std::string buf(std::istreambuf_iterator<char>(is), (std::istreambuf_iterator<char>()));

		vertex_shader = glCreateShader(GL_VERTEX_SHADER);
		const GLchar* buffer = buf.c_str();
		GLint size = buf.length();
		glShaderSource(vertex_shader, 1, &buffer, &size);
		glCompileShader(vertex_shader);
	}

	GLint compilationError(0);
	glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &compilationError);
	if (compilationError != GL_TRUE)
	{
		GLint errorSize(0);
		glGetShaderiv(vertex_shader, GL_INFO_LOG_LENGTH, &errorSize);
		char *error = new char[errorSize + 1];
		glGetShaderInfoLog(vertex_shader, errorSize, &errorSize, error);
		error[errorSize] = '\0';
		std::cout << "Error vertex: " << error << std::endl;
		delete[] error;
	}

	std::ifstream is2(fragmentShaderFilename);
	if (is2) {
//		std::string buf(std::istreambuf_iterator<char>(is2), {});
		std::string buf(std::istreambuf_iterator<char>(is2), (std::istreambuf_iterator<char>()));

		fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
		const GLchar* buffer = buf.c_str();
		GLint size = buf.length();
		glShaderSource(fragment_shader, 1, &buffer, &size);
		glCompileShader(fragment_shader);
	}

	glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &compilationError);
	if (compilationError != GL_TRUE)
	{
		GLint errorSize(0);
		glGetShaderiv(fragment_shader, GL_INFO_LOG_LENGTH, &errorSize);
		char *error = new char[errorSize + 1];
		glGetShaderInfoLog(fragment_shader, errorSize, &errorSize, error);
		error[errorSize] = '\0';
		std::cout << "Error fragment: " << error << std::endl;
		delete[] error;
	}

	glAttachShader(m_program, vertex_shader);
	glAttachShader(m_program, fragment_shader);

	glLinkProgram(m_program);

	GLint isLinked = 0;
	glGetProgramiv(m_program, GL_LINK_STATUS, &isLinked);
	if (isLinked == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetProgramiv(m_program, GL_INFO_LOG_LENGTH, &maxLength);
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(m_program, maxLength, &maxLength, &infoLog[0]);

		for (auto c : infoLog) {
			std::cout << c;
		}
		glDeleteProgram(m_program);
		return;
	}

	glDetachShader(m_program, vertex_shader);
	glDetachShader(m_program, fragment_shader);
}
