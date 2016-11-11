#include "particle.h"

Particle::Particle(MeshAsset* particle_asset, glm::vec3 position, glm::vec3 velocity): instance(particle_asset), position(position), velocity(velocity)
{
}

Particle::~Particle()
{
}
