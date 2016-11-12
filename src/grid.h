#pragma once

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "particle.h"

struct Cell {
	std::vector<Particle*> particles;
};

class Grid
{
public:
	Grid(glm::vec3 corner = glm::vec3(0, 0, 0), float scale = 0.1);
	~Grid();

	bool insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end);
	void computeNeighbours();

private:
	glm::vec3 m_corner;
	float m_scale, m_invScale;

	std::vector<std::vector<std::vector<Cell> > > m_cells;

	int m_sizeX, m_sizeY, m_sizeZ;
};
