#pragma once

#include <set>
#include <memory>

#include "mesh.h"

struct Particle
{

	Particle(MeshAsset* particle_asset,const float &size, const float &rest_density, bool meshParticle = false, glm::vec3 position = glm::vec3(0,0,0), glm::vec3 velocity = glm::vec3(0,0,0));
	~Particle();

	MeshInstance instance;

	std::set<Particle*> neighbours;

	float size;
	float mass;
	float volume;
	glm::vec3 position;
	glm::vec3 velocity;
	float density;
	float pressure;

	glm::vec3 localForce;

	bool isMeshParticle;
};

