#include "particle.h"

Particle::Particle(MeshAsset* particle_asset, int index, const float& size, const float& rest_density, bool meshParticle, glm::vec3 position, glm::vec3 velocity) : instance(particle_asset, position, glm::vec3(0,0,0), size / 2),
id(index),
size(size),
position(position),
velocity(velocity),
density(rest_density),
mass(1.33 * 3.14159 * std::pow(0.5f * size, 3) * rest_density),
isMeshParticle(meshParticle)
{
}

Particle::~Particle()
{
}
