#include "scene.h"


Scene::Scene(const std::string& initial_configuration) :  m_meshAsset(std::string("../data/meshes/plane.obj"), false), 
m_surface(&m_meshAsset, glm::vec3(0,0,0), glm::vec3(0,0,0), 0.3f),
m_sampler(m_surface.getFaces()),
m_fluid(initial_configuration)
{
	std::set<glm::vec3, CompareVec3> samples = m_sampler.sampleMeshUniformly(m_fluid.getSimulationScale() / 2.0f);

	std::vector<glm::vec3> positions(samples.begin(), samples.end());
	m_fluid.addBoundaryParticles(positions);
}

Scene::~Scene()
{

}

void Scene::update()
{

	//if (m_fluid.getTotalSimulationTime() < 2.0 && m_fluid.getTotalSimulationTime() - insertionTime > 0.05) {
	//	insertionTime = m_fluid.getTotalSimulationTime();
	//	for (int x = -2; x < 2; x++) {
	//		for (int y = -2; y < 2; y++) {
	//			glm::vec3 pos(0, 0.5, 1.5);
	//			pos += glm::vec3(0, x * (1.0 - (float)(std::rand() % 100) / 100.0f) * 0.1, y *(1.0 - (float)(std::rand() % 100) / 100.0f)* 0.1);
	//			m_fluid.addParticle(pos, glm::vec3(1.0f, 0, -1));
	//		}
	//	}
	//}



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
