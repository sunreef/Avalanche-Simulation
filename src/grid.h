#pragma once

#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "particle.h"


enum GType {
  GBoundary = 0,
  GFluid = -1,
  GOverlap = 1
};

struct Cell {

//  Cell(){density = 1.255f, rsnow = 0.f;}// set density to air density

	std::vector<Particle*> fluidParticles;
	std::vector<Particle*> boundaryParticles;

  // for snow smoke
  glm::vec3 velocity, velocity1, force, vorticity;
  glm::vec3 tmpVelocity;
  float divergence, density, pressure, tmpDensity;//rsnow, 
  glm::vec3 center, corner;
  GType type; // 0 for boundary, -1 for fluid snow, 1 for overlapped region

  void print() {
    printf("velocity %.6f %.6f %.6f, force %.6f %.6f %.6f, pressure %.6f, density %.6f, tmpDensity, %.6f, type %d\n",
            velocity.x, velocity.y, velocity.z, force.x, force.y, force.z, pressure, density, tmpDensity, type);
  }

  glm::vec3 getVelocityAt(const glm::vec3 &p) {
    float dx = p.x, dy = p.y, dz = p.z;
    if (dx < 0.f or dx > 1.f or dy < 0.f or dy > 1.f or dz < 0.f or dz > 1.f)
      printf("ERROR:: in Cell::getVelocityAt() with dx = %.3f dy = %.3f dz = %.3f\n", dx, dy, dz);

    return glm::vec3(
              (1.f - dx) * velocity.x + dx * velocity1.x,
              (1.f - dy) * velocity.y + dy * velocity1.y,
              (1.f - dz) * velocity.z + dz * velocity1.z
            );
  }
  glm::vec3 getVelocityAtCenter() {
      return 0.5f * (velocity + velocity1);   
  }

};

class Grid
{
public:
	Grid(glm::vec3 corner = glm::vec3(0, 0, 0), float scale = 0.1);
	~Grid();
  void clear();

	bool insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end);
  std::vector<glm::vec3> getFluidSnowBoundary(int density);
  void computeFluidSmokeBoundary();
  void sphGridInteraction(const float timeStep);
  void computeSmokeForce();
  void addSmokeForce(float timeStep);
  void solvePressure(float timeStep);
  void solvePoisson();
  void correctVelocity(float timeStep);
  void updateDensityAndSpeed(float timeStep);
  std::vector<glm::vec4> adhereSnowSmoke(float timeStep);
  bool smaller(const glm::vec3 &p) {return p.x+1 < m_corner.x or p.y+1 < m_corner.y or p.z+1 < m_corner.z;}

	void computeNeighboursFluid();
	void computeNeighboursBoundary();
  int getSearchRadius() {return m_searchRadius;}

  glm::vec3 worldToGrid(const glm::vec3 &p) {return (p - m_corner) * m_invScale;}
  glm::vec3 worldToGridIdx(const glm::vec3 &p) {
    glm::vec3 ret = (p - m_corner) * m_invScale;
    return glm::vec3(std::floor(ret.x), std::floor(ret.y), std::floor(ret.z));
  }
  glm::vec3 gridToWorld(const glm::vec3 &p) {return m_corner + p * m_scale;}
  const Cell& at(const glm::vec3 &idx) {
    return m_cells[idx.x][idx.y][idx.z];
  }
  
private:
	glm::vec3 m_corner;
	float m_scale, m_invScale;

	std::vector<std::vector<std::vector<Cell> > > m_cells;

	int m_sizeX, m_sizeY, m_sizeZ;

  int m_searchRadius;

  float k_drag, k_lift, k_cg, k_vorticity, k_airDensity;
  float k_accuracy;
  int k_iteration;

};
