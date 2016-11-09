#include "mesh.h"

#define BUFFER_OFFSET(i) ((GLintptr)NULL + (i))

size_t MeshAsset::count_meshes = 0;

MeshAsset::MeshAsset(const std::string &filename)
{
	count_meshes++;
	m_vbo = count_meshes;
	m_vao = count_meshes;

	loadObj(filename);
	initVAO();
}

MeshAsset::~MeshAsset()
{
}

void MeshAsset::loadObj(const std::string & filename)
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
		
			std::string vertex;
			for (int v = 0; v < 3; v++) {
				int i, j, k;
				ss >> vertex;
				std::stringstream ss2(vertex);
				std::string index;
				std::getline(ss2, index, '/');
				i = std::stoi(index);

				std::getline(ss2, index, '/');
				k = std::stoi(index);

				std::getline(ss2, index, '/');
				j = std::stoi(index);
				Vertex point;

				//std::cout << i << " " << j << " " << k << std::endl;

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

	m_numberOfVertices = vertices.size();

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
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(0), 3 * vertices.size() * sizeof(float), &final_positions[0]);
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(3 * vertices.size() * sizeof(float)), 3 * vertices.size() * sizeof(float), &final_normals[0]);
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(6 * vertices.size() * sizeof(float)), 2 * vertices.size() * sizeof(float), &final_textures[0]);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void MeshAsset::initVAO()
{
	if (glIsVertexArray(m_vao)) {
		glDeleteVertexArrays(1, &m_vao);
	}
	glGenVertexArrays(1, &m_vao);

	glBindVertexArray(m_vao);

	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)(0));
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)(3 * m_numberOfVertices * sizeof(float)));
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)(6 * m_numberOfVertices * sizeof(float)));

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void MeshAsset::destroy()
{
	glDeleteBuffers(1, &m_vbo);
	glDeleteVertexArrays(1, &m_vao);
}

//////////////////////////////////////////////////////////

MeshInstance::MeshInstance(MeshAsset * asset, const glm::vec3 & pos, const glm::vec3 & angles, float scale)
{
	m_asset = asset;
	m_position = pos; 
	m_angles = angles;
	m_scale = scale;
}

void MeshInstance::draw(Program & prog, const glm::mat4 & view)
{
	updateModelMatrix();
	glm::mat4 meshModel = m_modelMatrix;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << meshModel[i][j] << " ";
		}
		std::cout << std::endl;
	}

	prog.loadModelViewMatrix(view * meshModel);

	glBindVertexArray(m_asset->m_vao);

	glDrawArrays(GL_TRIANGLES, 0, m_asset->m_numberOfVertices);

	glBindVertexArray(0);
}

void MeshInstance::updateModelMatrix()
{
	glm::mat4 meshModel;

	meshModel = glm::scale(meshModel, glm::vec3(m_scale));
	meshModel = glm::rotate(meshModel, m_angles.x, glm::vec3(1, 0, 0));
	meshModel = glm::rotate(meshModel, m_angles.y, glm::vec3(0, 1, 0));
	meshModel = glm::rotate(meshModel, m_angles.z, glm::vec3(0, 0, 1));
	meshModel = glm::translate(meshModel, m_position);

	m_modelMatrix = meshModel;
}
