#include "scene.h"


Scene::Scene(const std::string& initial_configuration): m_fluid(initial_configuration)
{
}

Scene::~Scene()
{

}

void Scene::update()
{
	m_fluid.nextStep();
}

void Scene::draw(const Program & prog, const glm::mat4 & view)
{
	m_fluid.draw(prog, view);
}

float Scene::getTime()
{
	return m_fluid.getTotalSimulationTime();
}
