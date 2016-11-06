#include "mesh.h"

#define BUFFER_OFFSET(i) ((GLintptr)NULL + (i))

Mesh::Mesh(int index)
{
	m_vbo = index;
	m_vao = index;
}

Mesh::~Mesh()
{
}

void Mesh::loadObj(std::string & filename)
{
	std::ifstream obj_file(filename);
	std::string line;

	std::vector<float> positions;
	std::vector<float> normals;
	std::vector<float> textures;

	std::vector<Vertex> vertices;

	while (std::getline(obj_file, line)) {
		std::stringstream ss(line);

		std::string word;
		ss >> word;
		if (word.compare("v") == 0) {
			float x, y, z;
			ss >> x >> y >> z;
			positions.push_back(x);
			positions.push_back(y);
			positions.push_back(z);
		}
		if (word.compare("vn") == 0) {
			float x, y, z;
			ss >> x >> y >> z;
			normals.push_back(x);
			normals.push_back(y);
			normals.push_back(z);
		}
		if (word.compare("vt") == 0) {
			float x, y;
			ss >> x >> y;
			textures.push_back(x);
			textures.push_back(y);
		}
		if (word.compare("f") == 0) {
			int i, j, k;
			std::string vertex;
			for (int v = 0; v < 3; v++) {
				ss >> vertex;
				std::stringstream ss2(vertex, '/');
				ss2 >> i >> j >> k;
				Vertex point;

				i--;
				j--;
				k--;

				point.x = positions[3 * i];
				point.y = positions[3 * i + 1];
				point.z = positions[3 * i + 2];

				point.n_x = normals[3 * j];
				point.n_y = normals[3 * j + 1];
				point.n_z = normals[3 * j + 2];

				point.t_x = textures[2 * k];
				point.t_y = textures[2 * k + 1];

				vertices.push_back(point);
			}
		}
	}

	std::vector<float> final_positions(3 * vertices.size());
	std::vector<float> final_normals(3 * vertices.size());
	std::vector<float> final_textures(2 * vertices.size());

	for (int v = 0; v < vertices.size(); v++) {
		final_positions[3 * v] = vertices[v].x;
		final_positions[3 * v + 1] = vertices[v].y;
		final_positions[3 * v + 2] = vertices[v].z;

		final_normals[3 * v] = vertices[v].n_x;
		final_normals[3 * v + 1] = vertices[v].n_y;
		final_normals[3 * v + 2] = vertices[v].n_z;

		final_textures[2 * v] = vertices[v].t_x;
		final_textures[2 * v + 1] = vertices[v].t_y;
	}

	if (glIsBuffer(m_vbo)) {
		glDeleteBuffers(1, &m_vbo);
	}
	glGenBuffers(1, &m_vbo);

	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

	glBufferData(GL_ARRAY_BUFFER, 8 * vertices.size() * sizeof(float), 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, BUFFER_OFFSET(0), 3 * vertices.size() * sizeof(float), &final_positions[0]);
	glBufferSubData(GL_ARRAY_BUFFER, BUFFER_OFFSET(3 * vertices.size() * sizeof(float)), 3 * vertices.size() * sizeof(float), &final_normals[0]);
	glBufferSubData(GL_ARRAY_BUFFER, BUFFER_OFFSET(6 * vertices.size() * sizeof(float)), 2 * vertices.size() * sizeof(float), &final_textures[0]);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Mesh::initVAO()
{
	if (glIsVertexArray(m_vao)) {
		glDeleteVertexArrays(1, &m_vao);
	}
	glGenVertexArrays(1, &m_vao);

	glBindVertexArray(m_vao);

	// TODO: Add the buffer calls once the shaders have been written.

	glBindVertexArray(0);


}

void Mesh::destroy()
{
	glDeleteBuffers(1, &m_vbo);
	glDeleteVertexArrays(1, &m_vao);
}
