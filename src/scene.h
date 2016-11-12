#pragma once

#include <vector>

#include "grid.h"
#include "particle.h"

class Scene
{
public:
	Scene();
	~Scene();


private:
	MeshAsset m_particleAsset;
	std::vector<Particle*> m_particles;
	Grid m_grid;
};

