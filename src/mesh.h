#pragma once

#include <GL/glew.h>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

struct Vertex {
	int index;
	float x, y, z;
	float n_x, n_y, n_z;
	float t_x, t_y;
};

class Mesh
{
public:
	Mesh(int index);
	~Mesh();


	void draw();

	void loadObj(std::string &filename);
	void initVAO();

	void destroy();

private:
	GLuint m_vbo;
	GLuint m_vao;
};

