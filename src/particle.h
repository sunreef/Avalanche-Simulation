#pragma once

#include <set>
#include <memory>

#include "mesh.h"

struct Particle
{

	Particle(MeshAsset* particle_asset, glm::vec3 position = glm::vec3(0,0,0), glm::vec3 velocity = glm::vec3(0,0,0));
	~Particle();

	MeshInstance instance;

	std::set<Particle*> neighbours;

	glm::vec3 position;
	glm::vec3 velocity;
	double density;
	double pressure;

	glm::vec3 localForce;
};

