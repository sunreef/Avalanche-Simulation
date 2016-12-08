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

struct Face {
	std::vector<glm::vec3> vertices;
};

class MeshAsset
{
	friend class MeshInstance;


public:
	MeshAsset(const std::string& filename, bool sample);
	~MeshAsset();

	void destroy();
	const Vertex& getMeshVertex(int id, int v) const;
	int getMeshSize() const;

	std::vector<Face>& getFaces();

private:
	static size_t count_meshes;

	GLuint m_vbo;
	GLuint m_vao;

	int m_numberOfVertices;

	void loadObj(const std::string &filename, bool sample);
	void initVAO();

	glm::vec3 m_position;
	glm::vec3 m_angles;
	float m_scale;

	glm::mat4 m_modelMatrix;

	std::vector<Vertex> m_vertices;
	std::vector<Face> m_faces;
};

class MeshInstance {
public:
	MeshInstance(MeshAsset* asset, const glm::vec3 &pos = glm::vec3(0, 0, 0), const glm::vec3 &angles = glm::vec3(0, 0, 0), float scale = 1.0f);
	void draw(const Program &prog, const glm::mat4 &view);

	void setPosition(const glm::vec3& position);
	void setAngles(const glm::vec3& angles);
	void setScale(float scale);
	void setColor(const glm::vec4& color);

	glm::vec3 getMeshVertex(int id, int v) const;
	int getMeshSize() const;

	std::vector<Face> getFaces();

private:
	MeshAsset* m_asset;

	glm::vec3 m_position;
	glm::vec3 m_angles;
	float m_scale;

	glm::mat4 m_modelMatrix;

	glm::vec4 m_color;

	void updateModelMatrix();
};

