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
	std::vector<Particle*> m_particles;
	Grid m_grid;
};

