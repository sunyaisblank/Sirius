#version 330 core
in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D screenTexture;

// Lens Flare Parameters
// These should be uniforms, but hardcoding reasonable defaults if not set
uniform float threshold = 0.9;
uniform float intensity = 1.5;

void main() {
    vec4 baseColor = texture(screenTexture, TexCoord);
    
    // 1. ANAMORPHIC STREAK GENERATION
    // Sample horizontally to create the "JJ Abrams" style streak
    // commonly associated with bright accretion disks (Interstellar style)
    
    vec3 streak = vec3(0.0);
    float weightSum = 0.0;
    
    // Wide horizontal blur
    int samples = 12;
    float spread = 0.005; // Distance between samples
    
    for(int i = -samples; i <= samples; ++i) {
        // Skip center (added via baseColor)
        if(i == 0) continue;
        
        float w = 1.0 - abs(float(i)) / float(samples + 1);
        vec2 offset = vec2(float(i) * spread, 0.0);
        
        // Sample with clamp to edge
        vec2 uv = clamp(TexCoord + offset, 0.0, 1.0);
        vec3 sampleColor = texture(screenTexture, uv).rgb;
        
        // Threshold: only bright pixels contribute to flare
        vec3 brightPart = max(sampleColor - vec3(threshold), vec3(0.0));
        
        // Enhance brightness of the streak
        streak += brightPart * w;
        weightSum += w;
    }
    
    if (weightSum > 0.0) {
        streak /= weightSum;
    }
    
    // 2. COLORIZE STREAK
    // Blue/Cyan tint for anamorphic look
    vec3 streakColor = streak * vec3(0.4, 0.6, 1.0) * intensity * 4.0;
    
    // 3. GLOW / BLOOM (Simplified)
    // Small gaussian blur for general glow
    vec3 glow = vec3(0.0);
    float glowWeight = 0.0;
    for(int x = -2; x <= 2; ++x) {
        for(int y = -2; y <= 2; ++y) {
            vec2 off = vec2(x,y) * 0.003;
            vec3 s = texture(screenTexture, clamp(TexCoord + off, 0.0, 1.0)).rgb;
            float w = 1.0;
            vec3 b = max(s - vec3(threshold), vec3(0.0));
            glow += b * w;
            glowWeight += w;
        }
    }
    glow /= glowWeight;
    vec3 glowColor = glow * vec3(1.0, 0.9, 0.8) * intensity * 2.0;

    // Combine
    // Tone mapping is applied later or assumed linear? 
    // Usually we output linear and screen handles gamma.
    // Here we just additively blend.
    
    FragColor = vec4(baseColor.rgb + streakColor + glowColor, 1.0);
}
