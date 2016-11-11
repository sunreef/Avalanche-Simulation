#pragma once

#include "particle.h"

class Grid
{
public:
	Grid(glm::vec3 corner, double scale);
	~Grid();

private:
	glm::vec3 m_corner;
	double m_scale;
};

