#include "grid.h"


Grid::Grid(glm::vec3 corner, float scale) : m_corner(corner),
m_scale(scale),
m_invScale(1.0 / scale),
m_sizeX(1),
m_sizeY(1),
m_sizeZ(1)
{
	m_cells = std::vector<std::vector<std::vector<Cell> > >(1, std::vector<std::vector<Cell> >(1, std::vector<Cell>(1)));
}

Grid::~Grid()
{
}

bool Grid::insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end)
{
	for (auto it = begin; it != end; it++) {
		float x = (*it)->position.x;
		float y = (*it)->position.y;
		float z = (*it)->position.z;

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
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int x = 0; x < m_sizeX; x++) {
		for (int y = 0; y < m_sizeY; y++) {
			for (int z = 0; z < m_sizeZ; z++) {
				if (m_cells[x][y][z].particles.size() == 0) {
					continue;
				}
				std::set<Particle*> neighbours;
				int searchRadius = 1;

				while (neighbours.size() <= 1) {
					neighbours.clear();
					int xMin = std::max(0, x - searchRadius);
					int xMax = std::min(m_sizeX, x + searchRadius + 1);
					int yMin = std::max(0, y - searchRadius);
					int yMax = std::min(m_sizeY, y + searchRadius + 1);
					int zMin = std::max(0, z - searchRadius);
					int zMax = std::min(m_sizeZ, z + searchRadius + 1);

					for (int x2 = xMin; x2 < xMax; x2++) {
						for (int y2 = yMin; y2 < yMax; y2++) {
							for (int z2 = zMin; z2 < zMax; z2++) {
								neighbours.insert(std::begin(m_cells[x2][y2][z2].particles), std::end(m_cells[x2][y2][z2].particles));
							}
						}
					}
					searchRadius++;
				}

				for (Particle* p : m_cells[x][y][z].particles) {
					p->neighbours = neighbours;
				}
			}
		}
	}
}
