#include "fluid_engine.h"



FluidEngine::FluidEngine(const std::string& initial_configuration) : m_particleAsset(std::string("../data/meshes/sphere.obj"))
{
	std::ifstream init_file(initial_configuration);


	m_kernelSmoothingLength = 0.05f;
	m_gridResolution = 2 * m_kernelSmoothingLength;
	m_kernel = ExponentialKernel(m_kernelSmoothingLength);
	m_restDensity = 1000.0f;
	m_stiffness = 100.0f;
	m_viscosity = 0.01f;

	int n;
	init_file >> n;

	m_particles.reserve(n);

	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();

	for (int p = 0; p < n; p++) {
		float x, y, z;
		//init_file >> x >> y >> z;

		x = (float)(rand() % 1000) / 1000;
		y = (float)(rand() % 1000) / 1000;
		z = (float)(rand() % 1000) / 1000;

		minX = std::min(x, minX);
		minY = std::min(y, minY);
		minZ = std::min(z, minZ);
		m_particles.push_back(new Particle(&m_particleAsset, m_kernelSmoothingLength, m_restDensity, glm::vec3(x, y, z)));
	}

	minX -= EPSILON;
	minY -= EPSILON;
	minZ -= EPSILON;
	m_grid = Grid(glm::vec3(minX, minY, minZ), m_gridResolution);
}


FluidEngine::~FluidEngine()
{
	m_particleAsset.destroy();
	for (auto part : m_particles) {
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
	for (Particle* part : m_particles) {
		part->instance.draw(prog, view);
	}
}

glm::vec3 FluidEngine::kernelGradient(const glm::vec3 & xi, const glm::vec3& xj) const
{
	float norm = glm::distance(xi, xj);
	if (norm == 0) {
		return glm::vec3(0,0,0);
	}
	float kernelDeriv = m_kernel.evaluateDeriv(norm / m_kernelSmoothingLength);
	return (kernelDeriv / (norm * m_kernelSmoothingLength)) * (xi - xj);
}

void FluidEngine::buildGrid()
{
	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();

	for (int p = 0; p < m_particles.size(); p++) {
		glm::vec3 pos = m_particles[p]->position;
		minX = std::min(pos.x, minX);
		minY = std::min(pos.y, minY);
		minZ = std::min(pos.z, minZ);
	}

	minX -= EPSILON;
	minY -= EPSILON;
	minZ -= EPSILON;

	m_grid = Grid(glm::vec3(minX, minY, minZ), m_gridResolution);
	m_grid.insertParticles(m_particles.begin(), m_particles.end());
}

void FluidEngine::computeNeighbours()
{
	m_grid.computeNeighbours();
}

void FluidEngine::updateDensity()
{
	float invKernelSmoothingLength = 1.0f / m_kernelSmoothingLength;

#pragma omp parallel for
	for (int i = 0; i < m_particles.size(); i++) {
		Particle* p = m_particles[i];
		float density = 0.0f;
		for (Particle* neighbour : p->neighbours) {
			float norm = glm::distance(p->position, neighbour->position);
			density += neighbour->mass * m_kernel.evaluate(norm * invKernelSmoothingLength);
		}

		p->density = density;
		p->pressure = m_stiffness * (std::pow(p->density / m_restDensity, 7) - 1.0f);
	}
}

void FluidEngine::computeForce()
{
#pragma omp parallel for
	for (int i = 0; i < m_particles.size(); i++) {
		Particle* p = m_particles[i];
		glm::vec3 pressureForce(0, 0, 0);
		glm::vec3 viscosityForce(0, 0, 0);
		glm::vec3 otherForce(0, 0, 0);


		for (Particle* neighbour : p->neighbours) {
			glm::vec3 dir = p->position - neighbour->position;
			if (glm::length2(dir) < EPSILON) {
				continue;
			}
			glm::vec3 kernelGrad = kernelGradient(p->position, neighbour->position);


			pressureForce += neighbour->mass * (p->pressure / (p->density * p->density) + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;

			viscosityForce += (neighbour->mass * glm::dot(dir, kernelGrad) / (neighbour->density * glm::dot(dir, dir))) * (p->velocity - neighbour->velocity);
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

#pragma omp parallel for
	for (int i = 0; i < m_particles.size(); i++) {
		Particle* p = m_particles[i];
		p->velocity += m_timeStep * p->localForce;
		float norm = glm::length(p->velocity);
		newTimeStep = std::min(newTimeStep, 0.4f * m_kernelSmoothingLength / norm);
		p->position += m_timeStep * p->velocity;

		if (p->position.x < 0) {
			p->position.x = 0;
			p->velocity.x *= -1.0f;
		}
		if (p->position.y < 0) {
			p->position.y = 0;
			p->velocity.y *= -1;
		}
		if (p->position.z < 0) {
			p->position.z = 0;
			p->velocity.z *= -1;
		}

		if (p->position.x > 3) {
			p->position.x = 3;
			p->velocity.x *= -1;
		}
		if (p->position.y > 3) {
			p->position.y = 3;
			p->velocity.y *= -1;
		}
		if (p->position.z > 3) {
			p->position.z = 3;
			p->velocity.z *= -1;
		}
		p->instance.setPosition(p->position);
		p->instance.setColor(glm::vec3(norm/5.0, norm/5.0, 1.0));
	}
	m_timeStep = newTimeStep;
}
