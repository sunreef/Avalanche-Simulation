#pragma once

#include <algorithm>
#include <omp.h>

#include "particle.h"

struct Cell {
	std::vector<Particle*> particles;
};

class Grid
{
public:
	Grid(glm::vec3 corner, double scale);
	~Grid();

	bool insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end);

	void computeNeighbours();

private:
	glm::vec3 m_corner;
	double m_scale, m_invScale;

	std::vector<std::vector<std::vector<Cell> > > m_cells;
	
	int m_sizeX, m_sizeY, m_sizeZ;
};
