#pragma once

#include <cmath>
#include <algorithm>

#include <glm/gtx/norm.hpp>

#include "particle.h"
#include "grid.h"
#include "kernel.h"

#define EPSILON 0.000001f


class FluidEngine
{
public:
	FluidEngine(const std::string& initial_configuration);
	~FluidEngine();

	void addParticle(glm::vec3& position, glm::vec3& velocity);

	void basicSphNextStep();
	void pcisphNextStep();
	void draw(const Program &prog, const glm::mat4& view);
	float getTotalSimulationTime();

protected:
	MeshAsset m_particleAsset;
	MeshAsset m_surfaceAsset;
	MeshInstance* m_surface;

	std::vector<Particle*> m_fluidParticles;
	std::vector<Particle*> m_meshParticles;

	Grid m_grid;

	CubicKernel m_kernel;

	float m_totalTime = 0.0f;
	float m_kernelSmoothingLength = 0.1f;
	float m_gridResolution = 0.1f;
	float m_restDensity = 1.0f;
	float m_stiffness = 10000.0f;
	float m_viscosity = 0.000001f;
	float m_timeStep = 0.001f;

	glm::vec3 kernelGradient(const glm::vec3& xi, const glm::vec3& xj) const;

	void initializeEngine();
	void buildGrid();
	void computeNeighbours();

	void updateDensity();
	void computeViscosityForce();
	void computeExternalForce();

	void basicSphPressureForce();

	void pcisphInitializePressure();
	void pcisphPressureForce();


	void updatePositionAndSpeed();
	void computeMeshParticleMass();
};


class IISPH : public FluidEngine {
public:
	IISPH(const std::string& initial_configuration);

	void nextStep();
};

