#version 330

in vec4 frag_normal;

out vec4 frag_color;

uniform vec3 color;

void main() {
	float cosine = 0.2 + 5.0 * max(0,dot(frag_normal, vec4(0.0,0.0,1.0,0.0)));
	frag_color = vec4(color,1.0) * vec4(cosine, cosine, cosine, 1.0);
}