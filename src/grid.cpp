#include "grid.h"

float lerp(float a, float b, float ratio) {
  return a * (1.f - ratio) + b * ratio;
}

Grid::Grid(glm::vec3 corner, float scale) : m_corner(corner),
m_scale(scale),
m_invScale(1.0f / scale),
m_sizeX(1),
m_sizeY(1),
m_sizeZ(1)
{
	m_cells = std::vector<std::vector<std::vector<Cell> > >(1, std::vector<std::vector<Cell> >(1, std::vector<Cell>(1)));
  m_searchRadius = 1;
  k_accuracy = 1e-5 * pow(1, 2);
  k_drag = 0.81f, k_lift = 0.2f;
  k_cg = 0.2f;
  k_vorticity = 0.2f;
  k_iteration = 50;
  k_airDensity = 1.255f;
}

Grid::~Grid()
{
}

void Grid::clear() {
  for(int x = 0; x < m_sizeX; x++)
    for (int y = 0; y < m_sizeY; y++) 
      for (int z = 0; z < m_sizeZ; z++)
        m_cells[x][y][z].fluidParticles.clear(),
        m_cells[x][y][z].boundaryParticles.clear();
}

bool Grid::insertParticles(const std::vector<Particle*>::iterator & begin, const std::vector<Particle*>::iterator & end)
{
	for (auto it = begin; it != end; it++) {
		float x = (*it)->position.x;
		float y = (*it)->position.y;
		float z = (*it)->position.z;

		//std::cout << x << " " << y << " " << z << std::endl;
		//std::cout << m_cells.size() << std::endl;

		while (x >= m_corner.x + m_sizeX * m_scale - 0.1f) {
			m_cells.push_back(std::vector<std::vector<Cell> >(m_sizeY, std::vector<Cell>(m_sizeZ)));
			m_sizeX++;
		}
		while (y >= m_corner.y + m_sizeY * m_scale - 0.1f) {
			for (auto& vecX : m_cells) {
				vecX.push_back(std::vector<Cell>(m_sizeZ));
			}
			m_sizeY++;
		}
		while (z >= m_corner.z + m_sizeZ * m_scale - 0.1f) {
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
		if (!(*it)->isMeshParticle)
			m_cells[x_idx][y_idx][z_idx].fluidParticles.push_back(*it);
		else
			m_cells[x_idx][y_idx][z_idx].boundaryParticles.push_back(*it);
	}

  for (int x = 0; x < m_sizeX; x++)
    for (int y = 0; y < m_sizeY; y++)
      for (int z = 0; z < m_sizeZ; z++)
        m_cells[x][y][z].center = gridToWorld(glm::vec3((double)x+0.5,(double)y+0.5,(double)z+0.5)),
        m_cells[x][y][z].corner = gridToWorld(glm::vec3(x, y, z));


  return true;
}

std::vector<glm::vec3> Grid::getFluidSnowBoundary(int density) {
  std::vector<glm::vec3> ret;
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        if (m_cells[x][y][z].fluidParticles.size() == 0) {
          continue;
        } else if (m_cells[x][y][z].fluidParticles.size() < density) {
          float posx = x * m_scale + m_corner.x;
          float posy = y * m_scale + m_corner.y;
          float posz = (z+1) * m_scale + m_corner.z;
          ret.push_back(glm::vec3(posx, posy, posz));       
        }
      }
    }
  }  

  return ret;
}

void Grid::computeNeighboursBoundary()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int x = 0; x < m_sizeX; x++) {
		for (int y = 0; y < m_sizeY; y++) {
			for (int z = 0; z < m_sizeZ; z++) {
				if (m_cells[x][y][z].boundaryParticles.size() == 0) {
					continue;
				}
				std::set<Particle*> neighbours;

				neighbours.clear();
				int xMin = std::max(0, x - m_searchRadius);
				int xMax = std::min(m_sizeX, x + m_searchRadius + 1);
				int yMin = std::max(0, y - m_searchRadius);
				int yMax = std::min(m_sizeY, y + m_searchRadius + 1);
				int zMin = std::max(0, z - m_searchRadius);
				int zMax = std::min(m_sizeZ, z + m_searchRadius + 1);

				for (int x2 = xMin; x2 < xMax; x2++) {
					for (int y2 = yMin; y2 < yMax; y2++) {
						for (int z2 = zMin; z2 < zMax; z2++) {
							neighbours.insert(std::begin(m_cells[x2][y2][z2].boundaryParticles), std::end(m_cells[x2][y2][z2].boundaryParticles));
						}
					}
				}

				for (Particle* p : m_cells[x][y][z].boundaryParticles) {
					p->neighbours = neighbours;
				}
			}
		}
	}
}

void Grid::sphGridInteraction(const float timeStep) {
//  printf("sphGridIntersection start\n");
  // add drag force and lift force and generated snow smoke
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        if (m_cells[x][y][z].type != GOverlap) continue;

        int xp1 = std::min(m_sizeX-1, x+1), xm1 = std::max(0, x-1);
        int yp1 = std::min(m_sizeY-1, y+1), ym1 = std::max(0, y-1);
        int zp1 = std::min(m_sizeZ-1, z+1), zm1 = std::max(0, z-1);
        glm::vec3 dx = (m_cells[xp1][y][z].getVelocityAtCenter() - m_cells[xm1][y][z].getVelocityAtCenter());
        glm::vec3 dy = (m_cells[x][yp1][z].getVelocityAtCenter() - m_cells[x][ym1][z].getVelocityAtCenter());
        glm::vec3 dz = (m_cells[x][y][zp1].getVelocityAtCenter() - m_cells[x][y][zm1].getVelocityAtCenter());

        m_cells[x][y][z].vorticity = glm::vec3(dz.y - dy.z, dx.z - dz.x, dy.x - dx.y);

        glm::vec3 fLift, fDrag;
        float smoke = 0.f;

        for (Particle* p : m_cells[x][y][z].fluidParticles) {
          if (!p->isOverlapped) continue;
          glm::vec3 dv = p->velocity - m_cells[x][y][z].getVelocityAt(worldToGrid(p->position)-glm::vec3(x,y,z));
          fDrag += k_drag * p->mass / p->size * glm::length(dv) * dv;
          fLift += k_lift * p->mass * dv * m_cells[x][y][z].vorticity;
          float dsmoke = k_cg * (std::pow(std::abs(dv.x), 1.7) + std::pow(std::abs(dv.y), 1.7) + std::pow(std::abs(dv.z), 1.7)) * p->surface * timeStep;
//          printf("%.3f %.3f %.3f %.3f %.3f %.3f\n", k_cg, std::pow(dv.x, 2), std::pow(dv.y, 2), std::pow(dv.z, 2), p->surface, timeStep);
          p->mass -= dsmoke;
          smoke += dsmoke;
        }
        m_cells[x][y][z].force = fLift + fDrag;
        m_cells[x][y][z].density += smoke/std::pow(m_scale, 3);
//        if (m_cells[x][y][z].type == GOverlap and m_cells[x][y][z].density > 1e-5)
//          printf("%d %d %d : %.6f (%.6f / %.6f)\n", x, y, z, m_cells[x][y][z].density, smoke, std::pow(m_scale,3));
      }
    }
  }
}

void Grid::computeSmokeForce() {
//  printf("compute Smoke force start\n");
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        float density = (k_airDensity + m_cells[x][y][z].density);
        float rsnow = m_cells[x][y][z].density / density;
        glm::vec3 fBuoyance = -(1.f - rsnow) * density * glm::vec3(0, 0, -9.8);
        int xp1 = std::min(m_sizeX-1, x+1), xm1 = std::max(0, x-1);
        int yp1 = std::min(m_sizeX-1, y+1), ym1 = std::max(0, y-1);
        int zp1 = std::min(m_sizeX-1, z+1), zm1 = std::max(0, z-1);
        double dx = glm::length(m_cells[xp1][y][z].vorticity - m_cells[xm1][y][z].vorticity);// * m_invScale,
        double dy = glm::length(m_cells[x][yp1][z].vorticity - m_cells[x][ym1][z].vorticity);// * m_invScale,
        double dz = glm::length(m_cells[x][y][zp1].vorticity - m_cells[x][y][zm1].vorticity);// * m_invScale;
        glm::vec3 N(dx, dy, dz);
        N = N / (glm::length(N) + 1e-20f);// to avoid 0 norm
        glm::vec3 fConf = k_vorticity * m_scale * glm::cross(N, m_cells[x][y][z].vorticity);
        m_cells[x][y][z].force = m_cells[x][y][z].force + fBuoyance + fConf;
      }
    }
  }
}

void Grid::addSmokeForce(float timeStep) {
//  printf("add smoke force start\n");
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        m_cells[x][y][z].velocity = m_cells[x][y][z].velocity + timeStep * m_cells[x][y][z].force;
      }
    }
  }
}

void Grid::solvePressure(float timeStep) {
//  printf("solve pressure start\n");
  // set boundary velocity to 0
  // copy boundary 
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = m_sizeZ-1; z >= 0; z--) {

        if (x == 0 or x == m_sizeX - 1 or y == 0 or y == m_sizeY - 1 or z == 0 or z == m_sizeZ - 1 or
            m_cells[x][y][z].type != GOverlap)
          m_cells[x][y][z].velocity = m_cells[x][y][z].velocity1 = glm::vec3(0);

        if (m_cells[x][y][z].type == GFluid or m_cells[x][y][z].type == GBoundary)
          m_cells[x][y][z].velocity = m_cells[x][y][z].velocity1 = glm::vec3(0.f), 
          m_cells[x][y][z].pressure = m_cells[x][y][std::min(z+1, m_sizeZ-1)].pressure;

        if (x == 0) m_cells[x][y][z].pressure = m_cells[x+1][y][z].pressure;
        if (y == 0) m_cells[x][y][z].pressure = m_cells[x][y+1][z].pressure;
        if (z == 0) m_cells[x][y][z].pressure = m_cells[x][y][z+1].pressure;
        if (x == m_sizeX - 1) m_cells[x][y][z].pressure = m_cells[x-1][y][z].pressure;
        if (y == m_sizeY - 1) m_cells[x][y][z].pressure = m_cells[x][y-1][z].pressure;
        if (z == m_sizeZ - 1) m_cells[x][y][z].pressure = m_cells[x][y][z-1].pressure;
      }
    }
  }

  // calculate divergence
  float invScaleDt = m_invScale / timeStep;
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        glm::vec3 dv = m_cells[x][y][z].velocity1 - m_cells[x][y][z].velocity;
        m_cells[x][y][z].divergence = -(dv.x + dv.y + dv.z) * invScaleDt;
      }
    }
  }

  solvePoisson();
}

void Grid::solvePoisson() {
//  printf("solve poisson start\n");
  int iteration = k_iteration;
  float accuracy = k_accuracy;

  for (int it = 0; it < iteration; it++) {
    double residual = 0.f;
    float dx2 = m_scale*m_scale;
    int nGrid = 0;
    for (int z = 1; z < m_sizeZ - 1; z++) {
      for (int y = 1; y < m_sizeY - 1; y++) {
        for (int x = 1; x < m_sizeX - 1; x++) {
          if (m_cells[x][y][z].type == GFluid) continue;
          nGrid++;
          m_cells[x][y][z].pressure = 
            (dx2 * m_cells[x][y][z].divergence +
             m_cells[x-1][y][z].pressure + 
             m_cells[x+1][y][z].pressure + 
             m_cells[x][y+1][z].pressure + 
             m_cells[x][y-1][z].pressure + 
             m_cells[x][y][z+1].pressure + 
             m_cells[x][y][z-1].pressure) / 6.f;
          residual += std::pow(m_cells[x][y][z].pressure - m_cells[x][y][z].divergence, 2);
        }
      }
    }

    residual = sqrt(residual) / (nGrid);

    if(it == k_iteration - 1) 
      printf("Pressure solver: it=%d , res=%f \n", it, residual);
    if(residual < k_accuracy) {
      printf("Pressure solver: it=%d , res=%f, converged \n", it, residual);
      break;
    }
  }
}

void Grid::correctVelocity(float timeStep) {
//  printf("correct velocity start\n");
  for (int z = m_sizeZ - 2; z > 0; z--) {
    for (int y = m_sizeY - 2; y > 0; y--) {
      for (int x = m_sizeX - 2; x > 0; x--) {
        double dx = (m_cells[x][y][z].pressure - m_cells[x-1][y][z].pressure);
        double dy = (m_cells[x][y][z].pressure - m_cells[x][y-1][z].pressure);
        double dz = (m_cells[x][y][z].pressure - m_cells[x][y][z-1].pressure);
        m_cells[x][y][z].velocity = m_cells[x][y][z].velocity - timeStep * m_invScale * glm::vec3(dx, dy, dz);
        m_cells[x][y][z].velocity1 = glm::vec3(m_cells[x+1][y][z].velocity.x, m_cells[x][y+1][z].velocity.y, m_cells[x][y][z+1].velocity.z);
      }
    }
  }
}

void Grid::updateDensityAndSpeed(float timeStep) {
  
  printf("update density and spped start\n");
  // update x Velocity
  for (int x = 1; x < m_sizeX - 1; x++) {
    for (int y = 1; y < m_sizeY - 1; y++) {
      for (int z = 1; z < m_sizeZ - 1; z++) {
        glm::vec3 v = m_cells[x][y][z].getVelocityAtCenter();
        glm::vec3 p = m_cells[x][y][z].center - timeStep * v;
        glm::vec3 idx = worldToGridIdx(p);
        glm::vec3 c = m_cells[idx.x][idx.y][idx.z].center;
        glm::vec3 dx(1,0,0), dy(0,1,0), dz(0,0,1);
        m_cells[x][y][z].tmpVelocity.x = lerp(lerp(lerp(at(idx).velocity.x, at(idx+dx).velocity.x, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy).velocity.x, at(idx+dy+dx).velocity.x, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               lerp(lerp(at(idx+dz).velocity.x, at(idx+dx+dz).velocity.x, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy+dz).velocity.x, at(idx+dy+dx+dz).velocity.x, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               (p.z-c.z)*m_invScale);
      }
    }
  }

  // update y Velocity
  for (int x = 1; x < m_sizeX - 1; x++) {
    for (int y = 1; y < m_sizeY - 1; y++) {
      for (int z = 1; z < m_sizeZ - 1; z++) {
        glm::vec3 v = m_cells[x][y][z].getVelocityAtCenter();
        glm::vec3 p = m_cells[x][y][z].center - timeStep * v;
        glm::vec3 idx = worldToGridIdx(p);
        glm::vec3 c = m_cells[idx.x][idx.y][idx.z].center;
        glm::vec3 dx(1,0,0), dy(0,1,0), dz(0,0,1);
        m_cells[x][y][z].tmpVelocity.y = lerp(lerp(lerp(at(idx).velocity.y, at(idx+dx).velocity.y, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy).velocity.y, at(idx+dy+dx).velocity.y, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               lerp(lerp(at(idx+dz).velocity.y, at(idx+dx+dz).velocity.y, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy+dz).velocity.y, at(idx+dy+dx+dz).velocity.y, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               (p.z-c.z)*m_invScale);
      }
    }
  }


  // update z Velocity
  for (int x = 1; x < m_sizeX - 1; x++) {
    for (int y = 1; y < m_sizeY - 1; y++) {
      for (int z = 1; z < m_sizeZ - 1; z++) {
        glm::vec3 v = m_cells[x][y][z].getVelocityAtCenter();
        glm::vec3 p = m_cells[x][y][z].center - timeStep * v;
        glm::vec3 idx = worldToGridIdx(p);
        glm::vec3 c = m_cells[idx.x][idx.y][idx.z].center;
        glm::vec3 dx(1,0,0), dy(0,1,0), dz(0,0,1);
        m_cells[x][y][z].tmpVelocity.z = lerp(lerp(lerp(at(idx).velocity.z, at(idx+dx).velocity.z, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy).velocity.z, at(idx+dy+dx).velocity.z, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               lerp(lerp(at(idx+dz).velocity.z, at(idx+dx+dz).velocity.z, (p.x-c.x)*m_invScale), 
                                                    lerp(at(idx+dy+dz).velocity.z, at(idx+dy+dx+dz).velocity.z, (p.x-c.x)*m_invScale), 
                                                    (p.y-c.y)*m_invScale),
                                               (p.z-c.z)*m_invScale);
      }
    }
  }


  // updage density
  printf("update density %d %d %d\n", m_sizeX - 1, m_sizeY - 1, m_sizeZ - 1);
  for (int x = 1; x < m_sizeX - 1; x++) {
    for (int y = 1; y < m_sizeY - 1; y++) {
      for (int z = 1; z < m_sizeZ - 1; z++) {
        glm::vec3 v = m_cells[x][y][z].getVelocityAtCenter();
        glm::vec3 p = m_cells[x][y][z].center - timeStep * v;
        glm::vec3 idx = worldToGridIdx((p - glm::vec3(0.5,0.5,0.5)*m_scale));
        glm::vec3 c = m_cells[idx.x][idx.y][idx.z].center;
        glm::vec3 dx(1,0,0), dy(0,1,0), dz(0,0,1);
        m_cells[x][y][z].tmpDensity = lerp(lerp(lerp(at(idx).density, at(idx+dx).density, (p.x-c.x)*m_invScale), 
                                                lerp(at(idx+dy).density, at(idx+dy+dx).density, (p.x-c.x)*m_invScale), 
                                                (p.y-c.y)*m_invScale), 
                                           lerp(lerp(at(idx+dz).density, at(idx+dx+dz).density, (p.x-c.x)*m_invScale), 
                                                lerp(at(idx+dy+dz).density, at(idx+dy+dx+dz).density, (p.x-c.x)*m_invScale), 
                                                (p.y-c.y)*m_invScale), 
                                           (p.z-c.z)*m_invScale);
        glm::vec3 o = m_cells[x][y][z].center;
      }
    }
  }


  // update
  for (int x = 1; x < m_sizeX - 1; x++) {
    for (int y = 1; y < m_sizeY - 1; y++) {
      for (int z = 1; z < m_sizeZ - 1; z++) {
        if (m_cells[x][y][z].type != GOverlap) {
          m_cells[x][y][z].velocity = glm::vec3(0.f);
          m_cells[x][y][z].density = 0.f;
        } else {
          int xp1 = std::min(m_sizeX-1, x+1), xm1 = std::max(0, x-1);
          int yp1 = std::min(m_sizeX-1, y+1), ym1 = std::max(0, y-1);
          int zp1 = std::min(m_sizeX-1, z+1), zm1 = std::max(0, z-1);

          m_cells[x][y][z].velocity = m_cells[x][y][z].tmpVelocity;
          m_cells[x][y][z].velocity1 = glm::vec3(m_cells[xp1][y][z].tmpVelocity.x, m_cells[x][yp1][z].tmpVelocity.y, m_cells[x][y][zp1].tmpVelocity.z);
          m_cells[x][y][z].density = m_cells[x][y][z].tmpDensity;
        }
      }
    }
  }
}

std::vector<glm::vec4> Grid::adhereSnowSmoke(float timeStep) {

//  printf("adhere snow smoke start\n");
  std::vector<glm::vec4> smoke;
//  m_cells[26][22][20].print();
  for (int x = 0; x < m_sizeX; x++) {
    for (int y = 0; y < m_sizeY; y++) {
      for (int z = 0; z < m_sizeZ; z++) {
        if (m_cells[x][y][z].type != GOverlap or m_cells[x][y][z].density < 0.0001) continue;
//        printf("density %d %d %d : %.10f\n",x, y, z, m_cells[x][y][z].density);
//        if (m_cells[x][y][z].type == GOverlap)
        smoke.push_back(glm::vec4(m_cells[x][y][z].center.x, m_cells[x][y][z].center.y, m_cells[x][y][z].center.z,
                                    m_cells[x][y][z].density));
      }
    }
  }

  return smoke;
}

void Grid::computeFluidSmokeBoundary() {
  for (int x = 0; x < m_sizeX; x++)
    for (int y = 0; y < m_sizeY; y++)
      for (int z = 0; z < m_sizeZ; z++) {

        if (m_cells[x][y][z].fluidParticles.size() == 0 and m_cells[x][y][z].boundaryParticles.size() == 0)
          m_cells[x][y][z].type = GOverlap;

        bool flag = false;
        for (Particle* p : m_cells[x][y][z].fluidParticles) {
          if (not p->isOverlapped) {flag = true; break;}
        }

        m_cells[x][y][z].type = flag? GFluid : GOverlap;
      }

  for (int x = 1; x < m_sizeX - 1; x++)
    for (int y = 1; y < m_sizeY - 1; y++)
      for (int z = 1; z < m_sizeZ - 1; z++) {

        bool flag = false;

        for (int dx = -1; dx <= 1 and !flag; dx++)
          for (int dy = -1; dy <= 1 and !flag; dy++)
            for (int dz = -1; dz <= 1 and !flag; dz++) 
              flag |= m_cells[x+dx][y+dy][z+dz].type == GOverlap;     
        
        if (flag and m_cells[x][y][z].type == GFluid)
          m_cells[x][y][z].type = GBoundary;
      }
}

void Grid::computeNeighboursFluid()
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int x = 0; x < m_sizeX; x++) {
		for (int y = 0; y < m_sizeY; y++) {
			for (int z = 0; z < m_sizeZ; z++) {
				if (m_cells[x][y][z].fluidParticles.size() == 0) {
					continue;
				}
				std::set<Particle*> neighbours;

				neighbours.clear();
				int xMin = std::max(0, x - m_searchRadius);
				int xMax = std::min(m_sizeX, x + m_searchRadius + 1);
				int yMin = std::max(0, y - m_searchRadius);
				int yMax = std::min(m_sizeY, y + m_searchRadius + 1);
				int zMin = std::max(0, z - m_searchRadius);
				int zMax = std::min(m_sizeZ, z + m_searchRadius + 1);

				for (int x2 = xMin; x2 < xMax; x2++) {
					for (int y2 = yMin; y2 < yMax; y2++) {
						for (int z2 = zMin; z2 < zMax; z2++) {
							neighbours.insert(std::begin(m_cells[x2][y2][z2].fluidParticles), std::end(m_cells[x2][y2][z2].fluidParticles));
							neighbours.insert(std::begin(m_cells[x2][y2][z2].boundaryParticles), std::end(m_cells[x2][y2][z2].boundaryParticles));
						}
					}
				}

				for (Particle* p : m_cells[x][y][z].fluidParticles) {
					p->neighbours = neighbours;
				}
			}
		}
	}
}
