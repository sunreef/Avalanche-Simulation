#include "scene.h"


Scene::Scene(const std::string& initial_configuration) : m_fluid(initial_configuration)
{
}

Scene::~Scene()
{

}

void Scene::update()
{

	if (m_fluid.getTotalSimulationTime() < 2.0 && m_fluid.getTotalSimulationTime() - insertionTime > 0.01) {
		insertionTime = m_fluid.getTotalSimulationTime();
		for (int x = -2; x < 3; x++) {
			for (int y = -2; y < 3; y++) {
				glm::vec3 pos(0, 0.5, 3.0);
				pos += glm::vec3(0, x * 0.05, y * 0.05);
				m_fluid.addParticle(pos, glm::vec3(4, 0, 0));
			}
		}
	}



	m_fluid.basicSphNextStep();
	//m_fluid.pcisphNextStep();
	//m_fluid.nextStep();
}

void Scene::draw(const Program & prog, const glm::mat4 & view)
{
	m_fluid.draw(prog, view);
}

float Scene::getTime()
{
	return m_fluid.getTotalSimulationTime();
}
