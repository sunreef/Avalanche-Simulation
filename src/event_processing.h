#pragma once

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cmath>
#include <algorithm>

#define MY_PI 3.14159265359f



class EventProcessing {
public:
	int size_x;
	int size_y;


private:
	float speed = 12.0;
	float mouseSpeed = 0.1;
	float zoomSpeed = 5.0;

	float deltaTime;
	clock_t time;

	float x, y;
	float dx, dy;

	bool keys[349];
	bool mouse_buttons[8];

	glm::vec3 position;
	float horizontalAngle, verticalAngle;

	float fov;

public:
	EventProcessing() {
		deltaTime = 0; 
		time = clock();
		x = 0;
		y = 0;
		dx = 0;
		dy = 0;

		for (int i = 0; i < 349; i++) {
			keys[i] = false;
		}
		for (int i = 0; i < 8; i++) {
			mouse_buttons[i] = false;
		}

		position = glm::vec3(-2, -2, 2);
		horizontalAngle = MY_PI / 4;
		verticalAngle = - MY_PI / 9;

		fov = 45.0f;
	}

	void initHandler(GLFWwindow* window) {
		glfwGetWindowSize(window, &size_x, &size_y);
	}

	void mouseMotionCallback(GLFWwindow* window, double x_pos, double y_pos) {
		dx = x_pos - x;
		dy = y_pos - y;
		x = x_pos;
		y = y_pos;
	}

	void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
		if (action == GLFW_PRESS) {
			mouse_buttons[button] = true;
		}
		if (action == GLFW_RELEASE) {
			mouse_buttons[button] = false;
		}
	}

	void keyboardCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
		if (action == GLFW_PRESS) {
			keys[key] = true;
		}
		if (action == GLFW_RELEASE) {
			keys[key] = false;
		}
	}

	void scrollCallback(GLFWwindow* window, double x_offset, double y_offset) {
		fov -= zoomSpeed * y_offset;

		if (fov < 10) {
			fov = 10.0f;
		}
		if (fov > 170) {
			fov = 170.0f;
		}
	}

	bool processEvents(glm::mat4 &view, glm::mat4 &proj) {
		deltaTime = std::min(0.02f, (float)(clock() - time) / CLOCKS_PER_SEC);
		time = clock();

		if (mouse_buttons[GLFW_MOUSE_BUTTON_1]) {
			horizontalAngle -= mouseSpeed * std::min(deltaTime, 0.03f) * dx;
			verticalAngle += mouseSpeed *  std::min(deltaTime, 0.03f) * dy;
			dx = 0;
			dy = 0;
		}

		glm::vec3 direction = glm::vec3(
			cos(verticalAngle) * sin(horizontalAngle),
			cos(verticalAngle) * cos(horizontalAngle),
			sin(verticalAngle)
		);

		glm::vec3 left = glm::vec3(
			sin(horizontalAngle - 3.14f / 2.0f),
			cos(horizontalAngle - 3.14f / 2.0f),
			0
		);

		glm::vec3 up = glm::cross(direction, left);

		if (keys[GLFW_KEY_W]) {
			position += deltaTime * speed * direction;
		}
		if (keys[GLFW_KEY_S]) {
			position -= deltaTime * speed *direction;
		}
		if (keys[GLFW_KEY_A]) {
			position += deltaTime * speed * left;
		}
		if (keys[GLFW_KEY_D]) {
			position -= deltaTime * speed * left;
		}
		if (keys[GLFW_KEY_R]) {
			position += deltaTime * speed * up;
		}
		if (keys[GLFW_KEY_F]) {
			position -= deltaTime * speed * up;
		}
		if (keys[GLFW_KEY_ESCAPE]) {
			return false;
		}

		view = glm::lookAt(position, position + direction, up);
		proj = glm::perspective(glm::radians(fov), (float)size_x / size_y, 0.1f, 100.0f);
		return true;
	}
};
