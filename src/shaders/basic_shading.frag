#version 330

in vec4 frag_normal;

out vec4 frag_color;

void main() {
	float cosine = gl_FragCoord.w * 10;
	frag_color = vec4(cosine, cosine, cosine , 1.0);

}