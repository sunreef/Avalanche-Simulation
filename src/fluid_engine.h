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

	void nextStep();
	void draw(const Program &prog, const glm::mat4& view);

private:
	MeshAsset m_particleAsset;
  MeshAsset m_surfaceAsset;
  MeshInstance* m_surface;

	std::vector<Particle*> m_particles;
  std::vector<Particle*> m_surface_particles;
	Grid m_grid;

	ExponentialKernel m_kernel;

	float m_kernelSmoothingLength = 0.1f;
	float m_gridResolution = 0.1f;
	float m_restDensity = 1.0f;
	float m_stiffness = 10000.0f;
	float m_viscosity = 0.1f;
	float m_timeStep = 0.005f;

	glm::vec3 kernelGradient(const glm::vec3& xi, const glm::vec3& xj) const;

	void buildGrid();
	void computeNeighbours();
	void updateDensity();
	void computeForce();
	void updatePositionAndSpeed();

};

