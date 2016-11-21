#include "particle.h"

Particle::Particle(MeshAsset* particle_asset, const float& size, const float& rest_density, bool meshParticle, glm::vec3 position, glm::vec3 velocity) : instance(particle_asset, position, glm::vec3(0,0,0), size / 3),
size(size),
position(position),
velocity(velocity),
density(rest_density),
mass(std::pow(0.66 * size, 3) * rest_density),
isMeshParticle(meshParticle)
{
}

Particle::~Particle()
{
}
