
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>

#include <time.h>
#include <chrono>
#include <thread>
#include <cmath>

#include "mesh.h"
#include "program.h"
#include "callback_struct.h"

#define MY_PI 3.14159265359f

int size_x = 1920;
int size_y = 1080;

float azerty_rules = true;

CallbackStruct callbackStructure;

void display() {}

void mouseCallback(int button, int state, int x, int y) {
	callbackStructure.x = x;
	callbackStructure.y = y;

	callbackStructure.mouse_buttons[button] = (state == GLUT_DOWN);
}

void passiveMotionCallback(int x, int y) {
	callbackStructure.dx = 0;
	callbackStructure.dy = 0;

	callbackStructure.x = x;
	callbackStructure.y = y;
}

void motionCallback(int x, int y) {
	callbackStructure.dx = x - callbackStructure.x;
	callbackStructure.dy = y - callbackStructure.y;
	callbackStructure.x = x;
	callbackStructure.y = y;
}

void keyboardCallback(unsigned char key, int x, int y) {
	callbackStructure.ascii_keys[key] = true;
}

void keyboardUpCallback(unsigned char key, int x, int y) {
	callbackStructure.ascii_keys[key] = false;
}

void initOpenGLContext(int argc, char** argv) {


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(size_x, size_y);
	glutCreateWindow("Hello OpenGL!");

	glutDisplayFunc(display);
	glutMouseFunc(mouseCallback);
	glutMotionFunc(motionCallback);
	glutPassiveMotionFunc(passiveMotionCallback);
	glutKeyboardFunc(keyboardCallback);
	glutKeyboardUpFunc(keyboardUpCallback);

	glewInit();
}

int main(int argc, char** argv) {

	initOpenGLContext(argc, argv);

	Program prog("../src/shaders/basic_shading.vert", "../src/shaders/basic_shading.frag");

	Mesh m(std::string("../data/meshes/mesh_0.obj"));



	float horizontalAngle = MY_PI;
	float verticalAngle = 0.0f;
	float initialFOV = 45.0f;

	float speed = 10.0f;
	float mouseSpeed = 0.1f;

	clock_t time = clock();
	float deltaTime = 0;

	int frameTime = CLOCKS_PER_SEC / 60;
	glm::vec3 position = glm::vec3(7, 10, 3.0);
	glm::mat4 modelView;
	glm::mat4 proj;

	glEnable(GL_DEPTH_TEST);
	while (true) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glutMainLoopEvent();
		deltaTime = (float)(clock() - time) / CLOCKS_PER_SEC;
		time = clock();
		horizontalAngle -= mouseSpeed * deltaTime * callbackStructure.dx;
		callbackStructure.dx = 0;
		verticalAngle += mouseSpeed * deltaTime * callbackStructure.dy;
		callbackStructure.dy = 0;

		glm::vec3 direction(
			cos(verticalAngle) * sin(horizontalAngle),
			cos(verticalAngle) * cos(horizontalAngle),
			sin(verticalAngle)
		);

		glm::vec3 left = glm::vec3(
			sin(horizontalAngle - 3.14f / 2.0f),
			cos(horizontalAngle - 3.14f / 2.0f),
			0
		);

		glm::vec3 up = -glm::cross(left, direction);

		if (callbackStructure.ascii_keys[(azerty_rules) ? 'z' : 'w']) {
			position += deltaTime * speed * direction;
			std::cout << deltaTime * speed << std::endl;
		}
		if (callbackStructure.ascii_keys['s']) {
			position -= deltaTime * speed *direction;
		}
		if (callbackStructure.ascii_keys['d']) {
			position -= deltaTime * speed * left;
		}
		if (callbackStructure.ascii_keys[(azerty_rules) ? 'q' : 'a']) {
			position += deltaTime * speed * left;
		}
		if (callbackStructure.ascii_keys['r']) {
			position += deltaTime * speed * up;
		}
		if (callbackStructure.ascii_keys['f']) {
			position -= deltaTime * speed * up;
		}
		if (callbackStructure.ascii_keys[27]) {
			break;
		}

		modelView = glm::lookAt(position, position + direction, up);
		proj = glm::perspective(initialFOV, (float)size_x / size_y, 0.01f, 1000.0f);

		prog.useProgram();
		prog.loadModelViewMatrix(modelView);
		prog.loadProjMatrix(proj);

		m.draw();

		prog.stopUseProgram();
		glutSwapBuffers();
	}


	return 0;
}