#include "mesh.h"

#define BUFFER_OFFSET(i) ((GLintptr)NULL + (i))

size_t MeshAsset::count_meshes = 0;

MeshAsset::MeshAsset(const std::string &filename, bool sample)
{
	count_meshes++;
	m_vbo = count_meshes;
	m_vao = count_meshes;

	loadObj(filename, sample);
	initVAO();
}

MeshAsset::~MeshAsset()
{
}

void MeshAsset::loadObj(const std::string & filename, bool sample)
{
	std::ifstream obj_file(filename);
	std::string line;

	std::vector<float> positions;
	std::vector<float> normals;
	std::vector<float> textures;

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

				m_vertices.push_back(point);
			}
		}
	}

  if (sample) {
    return;
  }

	m_numberOfVertices = m_vertices.size();

	std::vector<float> final_positions(3 * m_vertices.size());
	std::vector<float> final_normals(3 * m_vertices.size());
	std::vector<float> final_textures(2 * m_vertices.size());

	for (int v = 0; v < m_vertices.size(); v++) {
		final_positions[3 * v] = m_vertices[v].x;
		final_positions[3 * v + 1] = m_vertices[v].y;
		final_positions[3 * v + 2] = m_vertices[v].z;

		final_normals[3 * v] = m_vertices[v].n_x;
		final_normals[3 * v + 1] = m_vertices[v].n_y;
		final_normals[3 * v + 2] = m_vertices[v].n_z;

		final_textures[2 * v] = m_vertices[v].t_x;
		final_textures[2 * v + 1] = m_vertices[v].t_y;
	}

	if (glIsBuffer(m_vbo)) {
		glDeleteBuffers(1, &m_vbo);
	}
	glGenBuffers(1, &m_vbo);

	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);

	glBufferData(GL_ARRAY_BUFFER, 8 * m_vertices.size() * sizeof(float), 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(0), 3 * m_vertices.size() * sizeof(float), &final_positions[0]);
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(3 * m_vertices.size() * sizeof(float)), 3 * m_vertices.size() * sizeof(float), &final_normals[0]);
	glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)(6 * m_vertices.size() * sizeof(float)), 2 * m_vertices.size() * sizeof(float), &final_textures[0]);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

  if (!sample) m_vertices.clear();
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

const Vertex& MeshAsset::getMeshVertex(int id, int v) const {
  if ( 3 * id + v >= m_vertices.size()) {
    fprintf(stderr, "[Error] MeshAsset::getMesh() input index %d*3+%d larger then max index %d\n", id, v, m_vertices.size());
    return Vertex();
  }
  return m_vertices[3 * id + v];
}

int MeshAsset::getMeshSize() const {
  return m_vertices.size() / 3;
}

//////////////////////////////////////////////////////////

MeshInstance::MeshInstance(MeshAsset * asset, const glm::vec3 & pos, const glm::vec3 & angles, float scale)
{
	m_asset = asset;
	m_position = pos;
	m_angles = angles;
	m_scale = scale;
	m_color = glm::vec3((float)(rand() % 1000) / 1000, (float)(rand() % 1000) / 1000, (float)(rand() % 1000) / 1000);
}

void MeshInstance::draw(const Program & prog, const glm::mat4 & view)
{
	updateModelMatrix();
	glm::mat4 meshModel = m_modelMatrix;

	prog.loadColorUniform(m_color);

	prog.loadModelViewMatrix(view * meshModel);

	glBindVertexArray(m_asset->m_vao);

	glDrawArrays(GL_TRIANGLES, 0, m_asset->m_numberOfVertices);

	glBindVertexArray(0);
}

glm::vec3 MeshInstance::getMeshVertex(int id, int v) const{
  const Vertex& vertex = m_asset->getMeshVertex(id, v);
  glm::vec3 pos(vertex.x, vertex.y, vertex.z);
  return pos * m_scale + m_position[v];
}

int MeshInstance::getMeshSize() const{
  return m_asset->getMeshSize();
}

void MeshInstance::setPosition(const glm::vec3 & position)
{
	m_position = position;
}

void MeshInstance::setAngles(const glm::vec3 & angles)
{
	m_angles = angles;
}

void MeshInstance::setScale(float scale)
{
	m_scale = scale;
}

void MeshInstance::setColor(const glm::vec3 & color)
{
	m_color = color;
}

void MeshInstance::updateModelMatrix()
{
	glm::mat4 meshModel;
	meshModel = glm::translate(meshModel, m_position);
	meshModel = glm::scale(meshModel, glm::vec3(m_scale));
	meshModel = glm::rotate(meshModel, glm::radians(m_angles.x), glm::vec3(1, 0, 0));
	meshModel = glm::rotate(meshModel, glm::radians(m_angles.y), glm::vec3(0, 1, 0));
	meshModel = glm::rotate(meshModel, glm::radians(m_angles.z), glm::vec3(0, 0, 1));


	m_modelMatrix = meshModel;
}
