#include "scene.h"



Scene::Scene(): m_grid(glm::vec3(0,0,0), 0.2), m_particleAsset(std::string("../data/meshes/sphere.obj"))
{
	for (int p = 0; p < 100000; p++) {
		float x = (float)(rand() % 5000) / 1000;
		float y = (float)(rand() % 5000) / 1000;
		float z = (float)(rand() % 5000) / 1000;
		m_particles.push_back(new Particle(&m_particleAsset, glm::vec3(x,y,z)));
	}

	m_grid.insertParticles(m_particles.begin(), m_particles.end());
	m_grid.computeNeighbours();
}


Scene::~Scene()
{
	m_particleAsset.destroy();
	for (auto part : m_particles) {
		delete part;
	}
}
