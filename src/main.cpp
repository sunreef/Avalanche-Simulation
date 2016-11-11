
#include <GL/glew.h>

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <time.h>
#include <chrono>
#include <thread>
#include <random>
#include <cmath>
#include <algorithm>

#include "mesh.h"
#include "program.h"
#include "event_processing.h"


EventProcessing eventHandler;

void mouseMotionCallback(GLFWwindow* window, double x_pos, double y_pos) {
	eventHandler.mouseMotionCallback(window, x_pos, y_pos);
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	eventHandler.mouseButtonCallback(window, button, action, mods);
}

void keyboardCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	eventHandler.keyboardCallback(window, key, scancode, action, mods);
}

void scrollCallback(GLFWwindow* window, double x_offset, double y_offset) {
	eventHandler.scrollCallback(window, x_offset, y_offset);
}

void errorCallback(int error, const char* description) {
    fprintf(stderr, "Error: %s\n", description);
}

int main(int argc, char** argv) {

	if (!glfwInit()) {
    fprintf(stderr, "[ERROR] Failed to init GLFW\n");
    exit;
  }
  glfwSetErrorCallback(errorCallback);

	glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
	glfwWindowHint(GLFW_SAMPLES, 16);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWmonitor* monitor = glfwGetPrimaryMonitor();
	const GLFWvidmode* mode = glfwGetVideoMode(monitor);
	GLFWwindow* window = glfwCreateWindow(mode->width, mode->height, "Avalanche simulation", monitor, NULL);
  if (window == NULL) {
    fprintf(stderr, "[ERROR] Failed to create GLFW window\n");
    exit;
  }

  int major, minor, rev;
  glfwGetVersion(&major, &minor, &rev);
  fprintf(stderr, "GLFW version : %d.%d.%d\n", major, minor, rev);

	glfwMakeContextCurrent(window);

	glfwSetKeyCallback(window, keyboardCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetCursorPosCallback(window, mouseMotionCallback);
	glfwSetScrollCallback(window, scrollCallback);

	eventHandler.initHandler(window);

	glewInit();
  fprintf(stderr, "OpenGL version : %s\n", glGetString(GL_VERSION));

	Program prog("../src/shaders/basic_shading.vert", "../src/shaders/basic_shading.frag");

	MeshAsset mesh_asset(std::string("../data/meshes/plane.obj"));
	MeshAsset particle_asset(std::string("../data/meshes/sphere.obj"));

	MeshInstance instance_of_mesh(&mesh_asset);
	std::vector<MeshInstance> particles;

	for (int p = 0; p < 10000; p++) {
		float x = (float)(rand() % 5000) / 1000;
		float y = (float)(rand() % 5000) / 1000;
		float z = (float)(rand() % 5000) / 1000;
		particles.push_back(MeshInstance(&particle_asset, glm::vec3(x, y, z), glm::vec3(0, 0, 0), 0.05));
	}

	glm::mat4 view;
	glm::mat4 proj /*= glm::perspective(initialFOV, (float)size_x / size_y, 0.01f, 1000.0f)*/;

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	while (true) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glfwPollEvents();
		bool keep_rendering = eventHandler.processEvents(view, proj);

		prog.useProgram();

		prog.loadProjMatrix(proj);
		instance_of_mesh.draw(prog, view);

		for (auto instance : particles) {
			instance.draw(prog, view);
		}
		prog.stopUseProgram();
		glfwSwapBuffers(window);

		if (!keep_rendering) {
			break;
		}
	}

	prog.destroy();
	mesh_asset.destroy();
	particle_asset.destroy();
	return 0;
}
