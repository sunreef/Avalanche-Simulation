#pragma once

#include <set>
#include <memory>

#include "mesh.h"

struct Particle
{

	Particle(MeshAsset* particle_asset, int index, const float &size, const float &rest_density, bool meshParticle = false, glm::vec3 position = glm::vec3(0,0,0), glm::vec3 velocity = glm::vec3(0,0,0));
	~Particle();

	MeshInstance instance;

	int id;

	std::set<Particle*> neighbours;

	float size;
  float surface;
	float mass;
	float volume;
	glm::vec3 position;
	glm::vec3 velocity;
	float density;
	float pressure;

	glm::vec3 pressureForce;
	glm::vec3 viscosityForce;
	glm::vec3 externalForce;


	glm::vec3 temp_position;
	glm::vec3 temp_velocity;

	bool isMeshParticle;
  bool isOverlapped;
};

