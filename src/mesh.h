#pragma once

#include <GL/glew.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include "program.h"

struct Vertex {
	int index;
	float x, y, z;
	float n_x, n_y, n_z;
	float t_x, t_y;
};

class MeshAsset
{
	friend class MeshInstance;


public:
	MeshAsset(std::string & filename);
	~MeshAsset();

	void destroy();

private:
	static size_t count_meshes;

	GLuint m_vbo;
	GLuint m_vao;

	int m_numberOfVertices;

	void loadObj(std::string &filename);
	void initVAO();

	glm::vec3 m_position;
	glm::vec3 m_angles;
	float m_scale;

	glm::mat4 m_modelMatrix;
};

class MeshInstance {
public:
	MeshInstance(MeshAsset* asset, const glm::vec3 &pos = glm::vec3(0, 0, 0), const glm::vec3 &angles = glm::vec3(0, 0, 0), float scale = 1.0f);
	void draw(Program &prog, const glm::mat4 &view);

	void setPosition(const glm::vec3& position);
	void setAngles(const glm::vec3& angles);
	void setScale(float scale);

private:
	MeshAsset* m_asset;

	glm::vec3 m_position;
	glm::vec3 m_angles;
	float m_scale;

	glm::mat4 m_modelMatrix;

	glm::vec3 m_color;

	void updateModelMatrix();
};

