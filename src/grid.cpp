#include "grid.h"


Grid::Grid(glm::vec3 corner, double scale) : m_corner(corner), m_scale(scale), m_invScale(1.0 / scale), m_sizeX(1), m_sizeY(1), m_sizeZ(1)
{
	m_cells = std::vector<std::vector<std::vector<Cell> > >(1, std::vector<std::vector<Cell> >(1, std::vector<Cell>(1)));
}

Grid::~Grid()
{
}

bool Grid::insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end)
{
	for (auto it = begin; it != end; it++) {
		double x = (*it)->position.x;
		double y = (*it)->position.y;
		double z = (*it)->position.z;

		while (x >= m_corner.x + m_sizeX * m_scale) {
			m_cells.push_back(std::vector<std::vector<Cell> >(m_sizeY, std::vector<Cell>(m_sizeZ)));
			m_sizeX++;
		}
		while (y >= m_corner.y + m_sizeY * m_scale) {
			for (auto& vecX : m_cells) {
				vecX.push_back(std::vector<Cell>(m_sizeZ));
			}
			m_sizeY++;
		}
		while (z >= m_corner.z + m_sizeZ * m_scale) {
			for (auto& vecX : m_cells) {
				for (auto& vecY : vecX) {
					vecY.push_back(Cell());
				}
			}
			m_sizeZ++;
		}

		int x_idx = std::floor((x - m_corner.x) * m_invScale);
		int y_idx = std::floor((y - m_corner.y) * m_invScale);
		int z_idx = std::floor((z - m_corner.z) * m_invScale);

		if (x_idx < 0 || y_idx < 0 || z_idx < 0) {
			std::cout << "Error: a particle has been provided that is located before the corner of the grid." << std::endl;
			return false;
		}
		m_cells[x_idx][y_idx][z_idx].particles.push_back(*it);
	}
	return true;
}

void Grid::computeNeighbours()
{
#pragma omp parallel for
	for (int x = 0; x < m_sizeX; x++) {
		for (int y = 0; y < m_sizeY; y++) {
			for (int z = 0; z < m_sizeZ; z++) {
				if (m_cells[x][y][z].particles.size() == 0) {
					continue;
				}
				std::set<Particle*> neighbours;

				int xMin = std::max(0, x - 1);
				int xMax = std::min(m_sizeX, x + 2);
				int yMin = std::max(0, y - 1);
				int yMax = std::min(m_sizeY, y + 2);
				int zMin = std::max(0, z - 1);
				int zMax = std::min(m_sizeZ, z + 2);

				for (int x2 = xMin; x2 < xMax; x2++) {
					for (int y2 = yMin; y2 < yMax; y2++) {
						for (int z2 = zMin; z2 < zMax; z2++) {
							neighbours.insert(std::begin(m_cells[x2][y2][z2].particles), std::end(m_cells[x2][y2][z2].particles));
						}
					}
				}

				for (Particle* p : m_cells[x][y][z].particles) {
					p->neighbours = neighbours;
				}
			}
		}
	}
}
