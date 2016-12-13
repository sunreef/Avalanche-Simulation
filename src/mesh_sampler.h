#pragma once

#include "mesh.h"
#include <unordered_set>
#include <set>
#include <algorithm>

struct HashVec3 {
	size_t operator()(const glm::vec3& vec) const{
		std::hash<float> hash_float;
		return hash_float(vec[0]) + (hash_float(vec[1]) << 10) + (hash_float(vec[2]) << 20);
	}
};

struct CompareVec3 {
	bool operator()(const glm::vec3& a, const glm::vec3& b) {
		return (a[0] < b[0]) || ((a[0] == b[0]) && (a[1] < b[1])) || ((a[0] == b[0]) && (a[1] == b[1]) && (a[2] < b[2]));
	}
};

class MeshSampler
{
public:
	MeshSampler(const std::vector<Face>& faces);
	~MeshSampler();

	std::set<glm::vec3, CompareVec3> sampleMeshUniformly(float samplingScale);

private:
	std::vector<Face> m_faces;
};

