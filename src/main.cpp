
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
#include "scene.h"


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
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWmonitor* monitor = glfwGetPrimaryMonitor();
	const GLFWvidmode* mode = glfwGetVideoMode(monitor);

	GLFWwindow* window = glfwCreateWindow(mode->width * 0.8, mode->height * 0.8, "Avalanche simulation", NULL, NULL);
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

	Scene scene("../data/particle_configurations/config_1.txt");

	Program prog("../src/shaders/basic_shading.vert", "../src/shaders/basic_shading.frag");

	MeshAsset mesh_asset("../data/meshes/plane.obj", true);
	MeshInstance mesh_instance(&mesh_asset);

	glm::mat4 view;
	glm::mat4 proj;

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	while (true) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glfwPollEvents();
		bool keep_rendering = eventHandler.processEvents(view, proj);

		prog.useProgram();

		prog.loadProjMatrix(proj);

		clock_t time = clock();
		scene.update();
		std::cout << "Updated physics of the scene in " << (float)(clock() - time) / CLOCKS_PER_SEC << " seconds.\n";
		time = clock();

		if (glfwGetWindowAttrib(window, GLFW_FOCUSED)) {
			scene.draw(prog, view);
		}

		std::cout << "Rendered scene in " << (float)(clock() - time) / CLOCKS_PER_SEC << " seconds.\n";

		//mesh_instance.draw(prog, view);

		std::cout << "Total time: " << scene.getTime() << " seconds." << std::endl;

		prog.stopUseProgram();
		glfwSwapBuffers(window);

		if (!keep_rendering) {
			break;
		}
	}
	prog.destroy();
	glfwTerminate();
	return 0;
}
