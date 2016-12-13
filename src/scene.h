#pragma once

#include <vector>
#include <random>

#include "fluid_engine.h"
#include "mesh_sampler.h"

class Scene
{
public:
	Scene(const std::string& initial_configuration);
	~Scene();

	void update();
	void draw(const Program &prog, const glm::mat4& view);
	float getTime();

private:
	FluidEngine m_fluid;
	MeshAsset m_meshAsset;
	MeshInstance m_surface;
	MeshSampler m_sampler;

	float insertionTime = 0.0f;
};

