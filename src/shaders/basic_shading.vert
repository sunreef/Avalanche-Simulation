#version 330

layout(location=0) in vec3 position;
layout(location=1) in vec3 normal;
layout(location=2) in vec2 texture;

uniform mat4 modelView;
uniform mat4 proj;

out vec4 frag_normal;

void main() {
	gl_Position = proj * modelView * vec4(position, 1.0);
	vec3 norm_normal = normalize(normal);
	frag_normal = modelView * vec4(norm_normal,0.0);
}