#pragma once

#include <GL/glew.h>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

struct Vertex {
	int index;
	float x, y, z;
	float n_x, n_y, n_z;
	float t_x, t_y;
};

class Mesh
{
public:
	Mesh(std::string & filename);
	~Mesh();

	void draw();
	void destroy();

private:
	static size_t count_meshes;

	GLuint m_vbo;
	GLuint m_vao;

	int m_numberOfVertices;

	void loadObj(std::string &filename);
	void initVAO();
};

