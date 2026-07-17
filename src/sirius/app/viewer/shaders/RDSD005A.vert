// =============================================================================
// RDSD005A.vert - Bloom Post-Process Vertex Shader
// Component ID: RDSD005A
// =============================================================================
// Simple passthrough vertex shader for fullscreen quad
// =============================================================================

#version 330 core

layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;

out vec2 TexCoord;

void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
    TexCoord = aTexCoord;
}
