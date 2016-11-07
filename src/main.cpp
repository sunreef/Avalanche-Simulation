
#include <GL/glew.h>
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>


#include <time.h>
#include <chrono>
#include <thread>

#include "mesh.h"
#include "program.h"

void display() {}

int main(int argc, char** argv) {

	int size_x = 600;
	int size_y = 400;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);

	glutInitWindowPosition(0,0);
	glutInitWindowSize(size_x, size_y);
	glutCreateWindow("Hello OpenGL!");

	glutDisplayFunc(display);

	glewInit();


	Program prog("../src/shaders/basic_shading.vert", "../src/shaders/basic_shading.frag");

	Mesh m(std::string("../data/meshes/mesh_0.obj"));

	glEnable(GL_DEPTH_TEST);
	while (true) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glm::mat4 modelView = glm::lookAt(glm::vec3(20.0, 10.0, 3.0), glm::vec3(4, -9, 1), glm::vec3(0, 0, 1.0));
		glm::mat4 proj = glm::perspective(glm::radians(70.0f), (float)size_x / (float)size_y, 0.01f, 1000.0f);

		prog.useProgram();
		prog.loadModelViewMatrix(modelView);
		prog.loadProjMatrix(proj);

		m.draw();
		std::cout << "Drawing" << std::endl;


		prog.stopUseProgram();

		glutSwapBuffers();
		std::chrono::milliseconds timespan(16);
		std::this_thread::sleep_for(timespan);
	}


	return 0;
}