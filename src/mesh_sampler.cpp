#include "mesh_sampler.h"



MeshSampler::MeshSampler(const std::vector<Face>& faces) : m_faces(faces)
{
}


MeshSampler::~MeshSampler()
{
}

std::set<glm::vec3, CompareVec3> MeshSampler::sampleMeshUniformly(float samplingScale)
{
	std::set<glm::vec3, CompareVec3> samples;

	for (Face f : m_faces) {
		glm::vec3 dir0 = f.vertices[1] - f.vertices[0];
		glm::vec3 dir1 = f.vertices[2] - f.vertices[1];
		glm::vec3 dir2 = f.vertices[0] - f.vertices[2];

		float norm0 = glm::length(dir0);
		float norm1 = glm::length(dir1);
		float norm2 = glm::length(dir2);

		dir0 /= norm0;
		dir1 /= norm1;
		dir2 /= norm2;

		float maxNorm = std::max(norm0, std::max(norm1, norm2));

		glm::vec3 p;
		glm::vec3 vec1, vec2, vec3;
		float l1, l2, l3;

		if (maxNorm == norm0) {
			p = f.vertices[0];
			vec1 = dir0;
			vec2 = -dir2;
			vec3 = -dir1;

			l1 = norm0;
			l2 = norm2;
			l3 = norm1;
		}
		else if (maxNorm == norm1) {
			p = f.vertices[1];
			vec1 = dir1;
			vec2 = -dir0;
			vec3 = -dir2;

			l1 = norm1;
			l2 = norm0;
			l3 = norm2;
		}
		else {
			p = f.vertices[2];
			vec1 = dir2;
			vec2 = -dir1;
			vec3 = -dir0;

			l1 = norm2;
			l2 = norm1;
			l3 = norm0;
		}

		float cosine12 = glm::dot(vec1, vec2);
		float factor = samplingScale / cosine12;

		samples.insert(p);
		int i = 1;

		while (i * factor <= l2) {
			glm::vec3 p1 = p + (i * samplingScale) * vec1;
			glm::vec3 p2 = p + (i * factor) * vec2;
			float norm = glm::distance(p1, p2);
			int round = (int)std::floor(norm / samplingScale);

			float dist = norm / (round + 1);
			glm::vec3 dir = (p2 - p1) / norm;
			for (int a = 1; a <= round; a++) {
				samples.insert(p1 + a * dist * dir);
			}
			i++;
		}
		if (i * factor > l1) {
			continue;
		}

		float offset = l2 * cosine12;
		glm::vec3 q = p + l2 * vec2;
		float cosine13 = glm::dot(vec1, vec3);

		while (i * samplingScale <= l1) {
			glm::vec3 p1 = p + (i * samplingScale) * vec1;
			glm::vec3 p2 = q + (i * samplingScale - offset) / cosine13 * vec3;
			float norm = glm::distance(p1, p2);
			int round = (int)std::floor(norm / samplingScale);

			float dist = norm / (round + 1);
			glm::vec3 dir = (p2 - p1) / norm;
			for (int a = 0; a <= round + 1; a++) {
				samples.insert(p1 + a * dist * dir);
			}
			i++;
		}
	}

	for (Face f : m_faces) {
		glm::vec3 dir0 = f.vertices[1] - f.vertices[0];
		glm::vec3 dir1 = f.vertices[2] - f.vertices[1];
		glm::vec3 dir2 = f.vertices[0] - f.vertices[2];

		float norm0 = glm::length(dir0);
		float norm1 = glm::length(dir1);
		float norm2 = glm::length(dir2);

		dir0 /= norm0;
		dir1 /= norm1;
		dir2 /= norm2;

		//for (float t = 0; t <= norm0; t += samplingScale) {
		//	samples.insert()
		//}
	}

	return samples;
}
