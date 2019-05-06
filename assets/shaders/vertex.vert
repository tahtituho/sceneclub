#version 330 core
layout (location = 0) in vec3 aPos;

uniform float time;
varying float[12] sines;
varying float[12] coses;

void main() {
    for(int i = 0; i < sines.length(); i++) {
        sines[i] = sin(time / (i + 1));
        coses[i] = cos(time / (i + 1)); 
    }
    gl_Position = vec4(aPos.xyz, 1.0);
}