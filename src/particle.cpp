#include "particle.h"

Particle::Particle(MeshAsset* particle_asset, int index, const float& size, const float& rest_density, bool meshParticle, glm::vec3 position, glm::vec3 velocity) : instance(particle_asset, position, glm::vec3(0,0,0), size * (meshParticle ? 0.5f : 0.5f)),
id(index),
size(size * 0.5f),
position(position),
velocity(velocity),
density(rest_density),
mass(std::pow(0.5f * size, 3) * rest_density),
isMeshParticle(meshParticle)
{
  surface = size * size * size * 0.5f * 3.14f;
}

Particle::~Particle()
{
}
