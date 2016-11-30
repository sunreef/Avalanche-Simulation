#include "fluid_engine.h"
#include <string>
#include <random>



FluidEngine::FluidEngine(const std::string& initial_configuration) : m_particleAsset(std::string("../data/meshes/sphere.obj"), false), m_surfaceAsset(std::string("../data/meshes/plane.obj"), true)
{

	//m_fluidParticles = std::vector<Particle*>();
	std::ifstream init_file(initial_configuration);


	m_kernelSmoothingLength = 0.05f;
	m_gridResolution = 2 * m_kernelSmoothingLength;
	m_kernel = CubicKernel(m_kernelSmoothingLength);
	m_restDensity = 1000.0f;
	m_stiffness = 10000.0f;
	m_viscosity = 0.0005f;

	int n;
	init_file >> n;

	m_fluidParticles.reserve(n);

	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();

	//for (int p = 0; p < 0.15 * n; p++) {
	//	float x, y, z;
	//	//init_file >> x >> y >> z;

	//	x = (float)(rand() % 1000) / 1000;
	//	y = (float)(rand() % 1000) / 1000;
	//	z = (float)(rand() % 500) / 1000 + 0.2;

	//	minX = std::min(x, minX);
	//	minY = std::min(y, minY);
	//	minZ = std::min(z, minZ);
	//	m_fluidParticles.push_back(new Particle(&m_particleAsset, m_fluidParticles.size(), m_kernelSmoothingLength, m_restDensity, false, glm::vec3(x, y, z)));
	//}

	/*float samplingScale = m_kernelSmoothingLength / 2;
	for (int x = 0; x < 1.0 / samplingScale; x++) {
		for (int y = 0; y < 0.8 / samplingScale; y++) {
			float x_pos = x * samplingScale;
			float y_pos = y * samplingScale;
			m_meshParticles.push_back(new Particle(&m_particleAsset, m_meshParticles.size(), m_kernelSmoothingLength, m_restDensity, true, glm::vec3(x_pos, y_pos, 0)));
			m_meshParticles.push_back(new Particle(&m_particleAsset, m_meshParticles.size(), m_kernelSmoothingLength, m_restDensity, true, glm::vec3(x_pos, y_pos, -samplingScale)));
		}
	}*/

	minX -= EPSILON;
	minY -= EPSILON;
	minZ -= EPSILON;
	m_grid = Grid(glm::vec3(minX, minY, minZ), m_gridResolution);

	/**
	  * sample the mesh surface with particles
	  */

	  m_surface = new MeshInstance(&m_surfaceAsset, glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), .25f);
	  //n = (int) (2.0/m_kernelSmoothingLength);
	  //n = 1000;
	  m_meshParticles.reserve(n);

	  // get radom generator;
	  std::default_random_engine generator;
	  std::uniform_real_distribution<float> sampler(0.0, 1.0);

	  for (int p = 0; p < n; p++) {
	  	float x, y, z;

	  	int id = (int)(m_surface->getMeshSize() * (sampler(generator)));
	  	float su1 = sqrtf(sampler(generator));
	  	float u = 1.f - su1, v = su1 * sampler(generator);

	  	glm::vec3 pos = u * m_surface->getMeshVertex(id, 0) +
	  		v * m_surface->getMeshVertex(id, 1) +
	  		(1.f - u - v) * m_surface->getMeshVertex(id, 2);

	  	m_meshParticles.push_back(new Particle(&m_particleAsset, n, m_kernelSmoothingLength, m_restDensity, true, pos, glm::vec3(0,0,0)));
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

void FluidEngine::addParticle(glm::vec3 & position, const glm::vec3 & velocity)
{
	int n = m_fluidParticles.size();
	m_fluidParticles.push_back(new Particle(&m_particleAsset, n, m_kernelSmoothingLength, m_restDensity, false, position, velocity));
}

void FluidEngine::basicSphNextStep()
{
	buildGrid();
	computeNeighbours();
	updateDensity();
	basicSphPressureForce();
	computeViscosityForce();
	computeExternalForce();
	updatePositionAndSpeed();
}

void FluidEngine::pcisphNextStep()
{
	buildGrid();
	computeNeighbours();
	computeViscosityForce();
	computeExternalForce();

	pcisphInitializePressure();
	pcisphPressureForce();

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

//	MeshAsset tmp(std::string("../data/meshes/plane.obj"), false);
//	MeshInstance test(&tmp, glm::vec3(0,0,20), glm::vec3(0,0,0), .25f);
//	test.setColor(glm::vec3(0, 1.0, 1.0));
//	test.draw(prog, view);
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
	for (auto p : m_fluidParticles) {
		p->pressure = 0;
	}
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
		p->pressure = std::max(0.0f,m_stiffness * (p->density / m_restDensity - 1.0f));
	}
}

void FluidEngine::basicSphPressureForce()
{
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		glm::vec3 pressureForce(0, 0, 0);
		float factor = p->pressure / (p->density * p->density);

		if (p->density > 0) {
			for (Particle* neighbour : p->neighbours) {
				if (neighbour->position == p->position) {
					continue;
				}
				glm::vec3 kernelGrad = kernelGradient(p->position, neighbour->position);
				if (neighbour->isMeshParticle) {
					pressureForce += m_restDensity * neighbour->volume *  (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;
				}
				else {
					pressureForce += neighbour->mass * (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;
				}

			}
		}

		pressureForce *= -p->mass;

		p->pressureForce = pressureForce;
	}
}

void FluidEngine::pcisphInitializePressure()
{
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		p->pressure = 0;
		p->pressureForce = glm::vec3(0, 0, 0);
	}
}

void FluidEngine::pcisphPressureForce()
{
	float error_threshold = 0.01 * m_restDensity;
	float min_iterations = 3;

	float invKernelSmoothingLength = 1.0f / m_kernelSmoothingLength;

	float density_error = 0.0f;
	int iterations = 0;
	float mass = m_fluidParticles[0]->mass;
	float beta_constant_factor = 2 * m_timeStep * m_timeStep * mass * mass / (m_restDensity * m_restDensity);


	float delta = 0.0f;
	glm::vec3 sum_kernel_grad(0, 0, 0);
	float sum_kernel_sq_norm = 0.0f;
	for (int x = -3; x <= 3; x++) {
		for (int y = -3; y <= 3; y++) {
			for (int z = -3; z <= 3; z++) {
				glm::vec3 false_pos(x * m_kernelSmoothingLength, y * m_kernelSmoothingLength, z * m_kernelSmoothingLength);
				glm::vec3 kernelGrad = kernelGradient(glm::vec3(0, 0, 0), false_pos);
				sum_kernel_grad += kernelGrad;
				sum_kernel_sq_norm += glm::dot(kernelGrad, kernelGrad);
			}
		}
	}
	delta = 1.0f / (beta_constant_factor * (glm::dot(sum_kernel_grad, sum_kernel_grad) + sum_kernel_sq_norm));
	std::cout << delta << std::endl;

	while (iterations < 10 && (density_error > error_threshold || iterations < min_iterations)) {
		density_error = 0.0f;
		std::cout << iterations << std::endl;
		//std::vector<glm::vec3> temp_positions(m_fluidParticles.size());
		//std::vector<glm::vec3> temp_velocities(m_fluidParticles.size());

#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			Particle* p = m_fluidParticles[i];
			p->temp_velocity = p->velocity + (m_timeStep / p->mass) * (p->viscosityForce + p->externalForce + p->pressureForce);
			p->temp_position = p->position + m_timeStep * p->temp_velocity;
			//std::cout << iterations << " " << glm::length(p->temp_velocity) << std::endl;
		}

#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			Particle* p = m_fluidParticles[i];

			float density = 0.0f;
			for (Particle* neighbour : p->neighbours) {
				float norm = glm::distance(p->temp_position, neighbour->temp_position);
				density += neighbour->mass * m_kernel.evaluate(norm * invKernelSmoothingLength);
			}
			p->density = density;
			//std::cout << density << std::endl;
			float density_variation = density - m_restDensity;
			if (density_variation > 900 || density_variation < -900) {
				//std::cout << "Bad density " << density - p->mass * m_kernel.evaluate(0) << " " << density_variation << " " << density_error << std::endl;
			}
			//std::cout << density_variation << std::endl;
			density_error = std::max(density_error, std::abs(density_variation));

			glm::vec3 sum_kernel_grad(0, 0, 0);
			float sum_kernel_sq_norm = 0.0f;
			for (Particle* neighbour : p->neighbours) {
				glm::vec3 kernelGrad = kernelGradient(p->temp_position, neighbour->temp_position);
				sum_kernel_grad += kernelGrad;
				sum_kernel_sq_norm += glm::dot(kernelGrad, kernelGrad);
			}
			//std::cout << delta << " " << 1.0f / (beta_constant_factor * (glm::dot(sum_kernel_grad, sum_kernel_grad) + sum_kernel_sq_norm)) << std::endl;
			//p->pressure += /*delta **/ density_variation / (beta_constant_factor * (glm::dot(sum_kernel_grad, sum_kernel_grad) + sum_kernel_sq_norm));
			p->pressure = m_stiffness * (density_variation / m_restDensity);
			//std::cout << p->pressure << std::endl;
		}

#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			Particle* p = m_fluidParticles[i];
			glm::vec3 pressureForce(0, 0, 0);
			float factor = 2.0f * p->pressure / (p->density * p->density);

			if (p->density > 0) {
				for (Particle* neighbour : p->neighbours) {
					glm::vec3 kernelGrad = kernelGradient(p->temp_position, neighbour->temp_position);
					/*std::cout <<iterations << " " << factor << " " << glm::length(kernelGrad) << std::endl;*/
					pressureForce += kernelGrad;
					/*	if (neighbour->isMeshParticle) {
							pressureForce += (factor) * kernelGrad;
						}
						else {
							pressureForce += (factor + neighbour->pressure / (neighbour->density * neighbour->density)) * kernelGrad;
						}*/

				}
			}

			pressureForce *= -p->mass * p->mass * factor;
			//std::cout << glm::length(pressureForce) << std::endl;
			p->pressureForce = pressureForce;
		}
		iterations++;
		std::cout << density_error << std::endl;
	}
}

void FluidEngine::computeViscosityForce()
{
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		glm::vec3 viscosityForce(0, 0, 0);

		if (p->density > 0) {
			for (Particle* neighbour : p->neighbours) {
				if (neighbour->position == p->position) {
					continue;
				}
				glm::vec3 dir = p->position - neighbour->position;
				glm::vec3 kernelGrad = kernelGradient(p->position, neighbour->position);

				if (neighbour->isMeshParticle) {
				}
				else {
					viscosityForce += (neighbour->mass * glm::dot(dir, kernelGrad) / (neighbour->density * glm::dot(dir, dir))) * (p->velocity - neighbour->velocity);
				}

			}
		}
		viscosityForce *= 2 * m_viscosity * p->mass;

		p->viscosityForce = viscosityForce;
	}
}

void FluidEngine::computeExternalForce()
{
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		glm::vec3 otherForce(0, 0, 0);

		otherForce += p->mass * glm::vec3(0, 0, -9.81);

		p->externalForce = otherForce;
	}
}

void FluidEngine::updatePositionAndSpeed()
{
	float newTimeStep = 0.01f;
	float reboundRatio = 0.8f;

#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		p->velocity += (m_timeStep / p->mass) * (p->pressureForce + p->viscosityForce + p->externalForce);
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
		if (p->position.z < 0) {
			p->position.z = 0;
			p->velocity.z *= -reboundRatio;
		}

		if (p->position.x > 1.01) {
			p->position.x = 1.01;
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
		part->temp_position = part->position;
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

IISPH::IISPH(const std::string & initial_configuration) : FluidEngine(initial_configuration) {
}


void IISPH::nextStep()
{
	float invSmoothingDistance = 1.0f / m_kernelSmoothingLength;
	float omega = 0.5f;

	std::vector<float> a_ii(m_fluidParticles.size());
	std::vector<float> rho_i_adv(m_fluidParticles.size());

	std::vector<glm::vec3> d_ii(m_fluidParticles.size());
	std::vector<glm::vec3> sum_dij_pjl(m_fluidParticles.size());


	buildGrid();
	computeNeighbours();


	//#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		p->density = 0;
		for (auto neighbour : p->neighbours) {
			p->density += neighbour->mass * m_kernel.evaluate(glm::distance(p->position, neighbour->position) * invSmoothingDistance);
		}
	}
	//#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		if (!p->isMeshParticle) {
			glm::vec3 viscosity_force(0, 0, 0);
			if (p->density > 0) {
				for (auto neighbour : p->neighbours) {
					if (neighbour->position == p->position) {
						continue;
					}
					glm::vec3 dir = p->position - neighbour->position;
					glm::vec3 kernelGrad = kernelGradient(p->position, neighbour->position);
					viscosity_force += (neighbour->mass * glm::dot(dir, kernelGrad) / (neighbour->density * glm::dot(dir, dir))) * (p->velocity - neighbour->velocity);
				}
			}
			viscosity_force *= 2 * m_viscosity * p->mass;
			glm::vec3 gravity_force(0, 0, -9.81 * p->mass);

			p->viscosityForce = viscosity_force;
			p->externalForce = gravity_force;
		}
	}
	//#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		if (!p->isMeshParticle)
			p->velocity += (m_timeStep / p->mass) * (p->viscosityForce + p->externalForce);
		std::cout << glm::length(p->velocity) << std::endl;

		glm::vec3 dii(0, 0, 0);
		for (Particle* neighbour : p->neighbours) {
			glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
			dii -= neighbour->mass * kernel_grad;
		}
		d_ii[p->id] = dii * (m_timeStep * m_timeStep / (p->density * p->density));
	}




#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		float rho = 0.0f;
		for (auto neighbour : p->neighbours) {
			glm::vec3 v_ij = p->velocity - neighbour->velocity;
			glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
			rho += neighbour->mass * glm::dot(v_ij, kernel_grad);
		}
		rho_i_adv[p->id] = m_timeStep * rho + p->density;
		p->pressure *= 0.5f;
	}
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		float aii = 0.0f;
		float factor = m_timeStep * m_timeStep * p->mass / (p->density * p->density);
		for (auto neighbour : p->neighbours) {
			glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
			glm::vec3 d_ji = factor * kernel_grad;

			aii += neighbour->mass * glm::dot(kernel_grad, d_ii[p->id] - d_ji);
		}
		a_ii[p->id] = aii;
	}


	int l = 0;
	float variation;
	while (true) {
		variation = 0.0f;
#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			Particle* p = m_fluidParticles[i];
			glm::vec3 sum(0, 0, 0);
			for (auto neighbour : p->neighbours) {
				glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
				sum -= (neighbour->mass * neighbour->pressure / (neighbour->density * neighbour->density)) * kernel_grad;
			}
			sum_dij_pjl[p->id] = m_timeStep * m_timeStep * sum;
		}

		std::vector<float> new_pressures(m_fluidParticles.size());
		//#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			Particle* p = m_fluidParticles[i];
			float sum = 0.0f;
			float factor = m_timeStep * m_timeStep * p->mass / (p->density * p->density);
			for (auto neighbour : p->neighbours) {
				glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
				glm::vec3 sum_dij = sum_dij_pjl[p->id];
				glm::vec3 djj_pjl = d_ii[neighbour->id] * neighbour->pressure;
				glm::vec3 d_ji = factor * kernel_grad;
				glm::vec3 sum_djk = sum_dij_pjl[neighbour->id] - d_ji * p->pressure;

				glm::vec3 total_vector = sum_dij - djj_pjl - sum_djk;
				sum += neighbour->mass * glm::dot(total_vector, kernel_grad);
				if (sum != sum) {
					std::cout << glm::length(sum_dij) << " " << glm::length(djj_pjl) << std::endl;
				}
			}
			float new_pressure = (1.0f - omega) * p->pressure + (omega / a_ii[p->id]) * (m_restDensity - rho_i_adv[p->id] - sum);
			std::cout << p->pressure << " " << 1.0f / a_ii[p->id] * (m_restDensity - rho_i_adv[p->id] - sum) << std::endl;
			variation = std::max(variation, std::abs((new_pressure - p->pressure) / new_pressure));
			new_pressures[p->id] = new_pressure;
			/*		if (new_pressure != new_pressure) {
						std::cout << p->pressure << " " << a_ii[i] << " " << sum << " " << rho_i_adv[i] << std::endl;
					}*/
		}
#pragma omp parallel for
		for (int i = 0; i < m_fluidParticles.size(); i++) {
			m_fluidParticles[i]->pressure = new_pressures[i];
		}
		std::cout << l << " " << variation << std::endl;
		if (variation < 0.01 && l >= 2) {
			break;
		}
		l++;
	}

#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		glm::vec3 pressure_force(0, 0, 0);

		float value = p->pressure / (p->density * p->density);
		for (auto neighbour : p->neighbours) {
			glm::vec3 kernel_grad = kernelGradient(p->position, neighbour->position);
			pressure_force += neighbour->mass * (value + neighbour->pressure / (neighbour->density * neighbour->density)) * kernel_grad;
		}
		pressure_force *= -p->mass;
		p->pressureForce = pressure_force;
	}
#pragma omp parallel for
	for (int i = 0; i < m_fluidParticles.size(); i++) {
		Particle* p = m_fluidParticles[i];
		float reboundRatio = 0.8f;
		if (!p->isMeshParticle) {
			p->velocity += (m_timeStep / p->mass) * p->pressureForce;
			p->position += m_timeStep * p->velocity;
			//std::cout << glm::length(p->position) << std::endl;

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
		}
	}
	m_totalTime += m_timeStep;
}

