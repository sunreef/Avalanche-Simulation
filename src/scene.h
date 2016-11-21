#pragma once

#include <vector>

#include "fluid_engine.h"

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
};

