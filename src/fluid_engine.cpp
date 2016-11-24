#include "fluid_engine.h"
#include <string>
#include <random>



FluidEngine::FluidEngine(const std::string& initial_configuration) : m_particleAsset(std::string("../data/meshes/sphere.obj"), false), m_surfaceAsset(std::string("../data/meshes/plane.obj"), true)
{
	std::ifstream init_file(initial_configuration);


	m_kernelSmoothingLength = 0.1f;
	m_gridResolution = 2 * m_kernelSmoothingLength;
	m_kernel = CubicKernel(m_kernelSmoothingLength);
	m_restDensity = 1000.0f;
	m_stiffness = 10000.0f;
	m_viscosity = 0.0001f;

	int n;
	init_file >> n;

	//m_fluidParticles.reserve(n);

	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();

	for (int p = 0; p < 0.1 * n; p++) {
		float x, y, z;
		//init_file >> x >> y >> z;

		x = (float)(rand() % 1000) / 1000;
		y = (float)(rand() % 1000) / 1000;
		z = (float)(rand() % 500) / 1000 + 0.2;

		minX = std::min(x, minX);
		minY = std::min(y, minY);
		minZ = std::min(z, minZ);
		m_fluidParticles.push_back(new Particle(&m_particleAsset, m_kernelSmoothingLength, m_restDensity, false, glm::vec3(x, y, z)));
	}

//	float samplingScale = m_kernelSmoothingLength / 2;
//	for (int x = 0; x < 1.0 / samplingScale; x++) {
//		for (int y = 0; y < 0.8 / samplingScale; y++) {
//			float x_pos = x * samplingScale;
//			float y_pos = y * samplingScale;
//			m_meshParticles.push_back(new Particle(&m_particleAsset, m_kernelSmoothingLength, m_restDensity, true, glm::vec3(x_pos, y_pos, 0)));
//			m_meshParticles.push_back(new Particle(&m_particleAsset, m_kernelSmoothingLength, m_restDensity, true, glm::vec3(x_pos, y_pos, -samplingScale)));
//		}
//	}

	minX -= EPSILON;
	minY -= EPSILON;
	minZ -= EPSILON;
	m_grid = Grid(glm::vec3(minX, minY, minZ), m_gridResolution);

/**
  * sample the mesh surface with particles
  */ 

  m_surface = new MeshInstance(&m_surfaceAsset, glm::vec3(0,0,0), glm::vec3(0, 0, 0), .25f);
  //n = (int) (2.0/m_kernelSmoothingLength);
  n = 1000;
  m_meshParticles.reserve(n);

  // get radom generator;
  std::default_random_engine generator;
  std::uniform_real_distribution<float> sampler(0.0,1.0);

  for (int p = 0; p < n; p++) {
    float x, y, z;

    int id = (int) (m_surface->getMeshSize() * (sampler(generator)));
    float su1 = sqrtf(sampler(generator));
    float u = 1.f - su1, v = su1 * sampler(generator);
    
    glm::vec3 pos = u * m_surface->getMeshVertex(id, 0) +
                    v * m_surface->getMeshVertex(id, 1) +
                    (1.f - u - v) * m_surface->getMeshVertex(id, 2);
    
    m_meshParticles.push_back(new Particle(&m_particleAsset, m_kernelSmoothingLength, m_restDensity, true, pos));
  }

	initializeEngine();
}


FluidEngine::~FluidEngine()
{
	m_particleAsset.destroy();
	for (auto part : m_fluidParticles) {
		delete part;
	}
}

void FluidEngine::nextStep()
{
	buildGrid();
	computeNeighbours();
	updateDensity();
	computeForce();
	updatePositionAndSpeed();
}

void FluidEngine::draw(const Program & prog, const glm::mat4 & view)
{
	for (Particle* part : m_fluidParticles) {
		part->instance.draw(prog, view);
	}
	for (Particle* part : m_meshParticles) {
		part->instance.setColor(glm::vec3(0.8, 0.8, 0));
		part->instance.draw(prog, view);
	}

//  MeshAsset tmp(std::string("../data/meshes/plane.obj"), false);
//  MeshInstance test(&tmp, glm::vec3(0,0,0), glm::vec3(0,0,0), .25f);
//  test.setColor(glm::vec3(0, 1.0, 1.0));
//  test.draw(prog, view);
}

float FluidEngine::getTotalSimulationTime()
{
	return m_totalTime;
}

glm::vec3 FluidEngine::kernelGradient(const glm::vec3 & xi, const glm::vec3& xj) const
{
	float norm = glm::distance(xi, xj);
	if (norm == 0) {
		return glm::vec3(0, 0, 0);
	}
	float kernelDeriv = m_kernel.evaluateDeriv(norm / m_kernelSmoothingLength);
	return (kernelDeriv / (norm * m_kernelSmoothingLength)) * (xi - xj);
}

void FluidEngine::initializeEngine()
{
	buildGrid();
	computeNeighbours();
	computeMeshParticleMass();
}

void FluidEngine::buildGrid()
{
	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();

	for (int p = 0; p < m_fluidParticles.size(); p++) {
		glm::vec3 pos = m_fluidParticles[p]->position;
		minX = std::min(pos.x, minX);
		minY = std::min(pos.y, minY);
		minZ = std::min(pos.z, minZ);
	}
	for (int p = 0; p < m_meshParticles.size(); p++) {
		glm::vec3 pos = m_meshParticles[p]->position;
		minX = std::min(pos.x, minX);
		minY = std::min(pos.y, minY);
		minZ = std::min(pos.z, minZ);
	}

	minX -= EPSILON;
	minY -= EPSILON;
	minZ -= EPSILON;

	m_grid = Grid(glm::vec3(minX, minY, minZ), m_gridResolution);
	m_grid.insertParticles(m_fluidParticles.begin(), m_fluidParticles.end());
	m_grid.insertParticles(m_meshParticles.begin(), m_meshParticles.end());

}

void FluidEngine::computeNeighbours()
{
	m_grid.computeNeighbours();
}

void FluidEngine::updateDensity()
{
	float invKernelSmoothingLength = 1.0f / m_kernelSmoothingLength;

#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		float density = 0.0f;
		for (Particle* neighbour : p->neighbours) {
			float norm = glm::distance(p->position, neighbour->position);
			density += neighbour->mass * m_kernel.evaluate(norm * invKernelSmoothingLength);
		}

		p->density = density;
		//p->pressure = m_stiffness * (std::pow(p->density / m_restDensity, 7) - 1.0f);
		p->pressure = m_stiffness * (p->density / m_restDensity - 1.0);
	}

#pragma omp parallel for
	for (int i = 0; i < m_meshParticles.size(); i++) {
		Particle* p = m_meshParticles[i];
		float density = m_restDensity;
		for (Particle* neighbour : p->neighbours) {
			if (!neighbour->isMeshParticle) {
				float norm = glm::distance(p->position, neighbour->position);
				density += neighbour->mass * m_kernel.evaluate(norm * invKernelSmoothingLength);
			}

		}

		p->density = density;
		//p->pressure = m_stiffness * (std::pow(p->density / m_restDensity, 7) - 1.0f);
		p->pressure = m_stiffness * (p->density / m_restDensity - 1.0);
	}
}

void FluidEngine::computeForce()
{
	//#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		glm::vec3 pressureForce(0, 0, 0);
		glm::vec3 viscosityForce(0, 0, 0);
		glm::vec3 otherForce(0, 0, 0);

		float factor = p->pressure / (p->density * p->density);

		if (p->density > 0) {
			for (Particle* neighbour : p->neighbours) {
				if (neighbour->position == p->position) {
					continue;
				}
				glm::vec3 dir = p->position - neighbour->position;
				/*		if (glm::length2(dir) < EPSILON) {
							continue;
						}*/
				glm::vec3 kernelGrad = kernelGradient(p->position, neighbour->position);

				if (neighbour->isMeshParticle) {
					pressureForce += m_restDensity * neighbour->volume *  (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;
					//float nnn = glm::length(m_restDensity * neighbour->volume * factor * kernelGrad);
					//if (nnn > 0)
					//std::cout <<nnn << std::endl;
				}
				else {
					/*			float aaa = glm::length(neighbour->mass * (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad);
								std::cout << aaa << std::endl;*/
					pressureForce += neighbour->mass * (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;
					viscosityForce += (neighbour->mass * glm::dot(dir, kernelGrad) / (neighbour->density * glm::dot(dir, dir))) * (p->velocity - neighbour->velocity);
				}

			}
		}

		pressureForce *= -1.0f;
		viscosityForce *= 2 * m_viscosity;
		otherForce += glm::vec3(0, 0, -9.81);

		p->localForce = pressureForce + viscosityForce + otherForce;
	}
}

void FluidEngine::updatePositionAndSpeed()
{
	float newTimeStep = 0.01f;
	float reboundRatio = 0.8f;

#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		p->velocity += m_timeStep * p->localForce;
		float norm = glm::length(p->velocity);
		newTimeStep = std::min(newTimeStep, 0.4f * m_kernelSmoothingLength / norm);
		p->position += m_timeStep * p->velocity;


		if (p->position.x < 0) {
			p->position.x = 0;
			p->velocity.x *= -reboundRatio;
		}
		if (p->position.y < 0) {
			p->position.y = 0;
			p->velocity.y *= -reboundRatio;
		}
		if (p->position.z < -p->position.x * 0.5) {
			p->position.z = -p->position.x * 0.5;
			p->velocity -= (glm::dot(glm::vec3(1, 0, 1), p->velocity)) *glm::vec3(1, 0, 1);
		}

		if (p->position.x > 2) {
			p->position.x = 2;
			p->velocity.x *= -reboundRatio;
		}
		if (p->position.y > 1.01) {
			p->position.y = 1.01;
			p->velocity.y *= -reboundRatio;
		}
		if (p->position.z > 5) {
			p->position.z = 10 - p->position.z;
			p->velocity.z *= -reboundRatio;
		}
		p->instance.setPosition(p->position);
		if (p->isMeshParticle) {
			p->instance.setColor(glm::vec3(0.9, 0.9, 0.0));
		}
		else {
			norm = std::min(norm, 7.0f);
			//float colorDensity = std::min(p->density / m_restDensity, 2.0f);
			p->instance.setColor(glm::vec3(norm / 7.0, norm / 7.0, 1.0f));
		}
	}
	m_totalTime += m_timeStep;
	m_timeStep = newTimeStep;
}

void FluidEngine::computeMeshParticleMass()
{
	float invKernelSmoothingLength = 1.0f / m_kernelSmoothingLength;

#pragma omp parallel for
	for (int p = 0; p < m_meshParticles.size(); p++) {
		Particle* part = m_meshParticles[p];
		float inv_volume = 0.0f;
		for (Particle* neighbour : part->neighbours) {
			if (neighbour->isMeshParticle) {
				glm::vec3 vec = part->position - neighbour->position;
				inv_volume += m_kernel.evaluate(glm::length(vec) * invKernelSmoothingLength);
				//std::cout << m_kernel.evaluate(glm::length(vec) * invKernelSmoothingLength) << std::endl;
			}

		}
		//std::cout << inv_volume << std::endl;
		part->volume = 1.0f / inv_volume;
		part->mass = m_restDensity * part->volume;
	}
}
