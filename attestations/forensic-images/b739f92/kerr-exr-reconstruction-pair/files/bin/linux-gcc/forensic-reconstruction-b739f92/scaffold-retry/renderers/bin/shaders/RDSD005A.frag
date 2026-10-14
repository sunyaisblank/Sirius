// =============================================================================
// RDSD005A.frag - Bloom Post-Process Fragment Shader
// Component ID: RDSD005A
// =============================================================================
// Implements cinematic bloom/glow effect for Interstellar-style accretion disk
// 
// Algorithm:
//   1. Extract bright pixels above threshold
//   2. Apply Gaussian blur (combined horizontal + vertical approximation)
//   3. Additively composite with original image
//
// REFERENCES:
//   - Kawase bloom (efficient multi-pass blur)
//   - DNGR Paper: Lens flare / PSF convolution concept
// =============================================================================

#version 330 core

in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D screenTexture;
uniform vec2 texelSize;           // 1.0 / resolution
uniform float bloomIntensity;     // Bloom strength (0.0 - 1.0)
uniform float bloomThreshold;     // Brightness threshold for extraction
uniform int bloomEnabled;         // Toggle (0 = off, 1 = on)

// Sample offsets for efficient blur (Kawase-style)
const int BLUR_SAMPLES = 9;
const vec2 offsets[BLUR_SAMPLES] = vec2[](
    vec2(-1.0, -1.0), vec2(0.0, -1.0), vec2(1.0, -1.0),
    vec2(-1.0,  0.0), vec2(0.0,  0.0), vec2(1.0,  0.0),
    vec2(-1.0,  1.0), vec2(0.0,  1.0), vec2(1.0,  1.0)
);

// Gaussian-like weights (sum close to 1.0)
const float weights[BLUR_SAMPLES] = float[](
    0.0625, 0.125, 0.0625,
    0.125,  0.25,  0.125,
    0.0625, 0.125, 0.0625
);

// Luminance for perceptual brightness
float luminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

// Extract bright parts above threshold
vec3 extractBright(vec3 color, float threshold) {
    float lum = luminance(color);
    float contrib = max(0.0, lum - threshold) / max(lum, 0.0001);
    return color * contrib;
}

// Multi-scale blur combining multiple sample radii
vec3 multiScaleBlur(sampler2D tex, vec2 uv, vec2 texel, float radius) {
    vec3 result = vec3(0.0);
    
    // First pass: tight blur
    for (int i = 0; i < BLUR_SAMPLES; i++) {
        vec2 offset = offsets[i] * texel * radius;
        result += texture(tex, uv + offset).rgb * weights[i];
    }
    
    // Second pass: wider blur for glow
    vec3 wide = vec3(0.0);
    float wideRadius = radius * 2.5;
    for (int i = 0; i < BLUR_SAMPLES; i++) {
        vec2 offset = offsets[i] * texel * wideRadius;
        wide += texture(tex, uv + offset).rgb * weights[i];
    }
    
    // Third pass: even wider for atmospheric effect
    vec3 veryWide = vec3(0.0);
    float veryWideRadius = radius * 5.0;
    for (int i = 0; i < BLUR_SAMPLES; i++) {
        vec2 offset = offsets[i] * texel * veryWideRadius;
        veryWide += texture(tex, uv + offset).rgb * weights[i];
    }
    
    // Combine with falloff
    return result * 0.5 + wide * 0.35 + veryWide * 0.15;
}

void main() {
    vec3 originalColor = texture(screenTexture, TexCoord).rgb;
    
    if (bloomEnabled == 0) {
        FragColor = vec4(originalColor, 1.0);
        return;
    }
    
    // Extract bright pixels from scene
    vec3 brightColor = extractBright(originalColor, bloomThreshold);
    
    // Apply multi-scale blur for soft glow
    vec3 blurredBright = multiScaleBlur(screenTexture, TexCoord, texelSize, 4.0);
    
    // Extract brightness from blurred result too
    blurredBright = extractBright(blurredBright, bloomThreshold * 0.5);
    
    // Additive blend: original + bloom glow
    vec3 finalColor = originalColor + blurredBright * bloomIntensity;
    
    // Soft tone mapping to prevent clipping
    finalColor = finalColor / (finalColor + vec3(1.0));
    
    // Re-apply some contrast
    finalColor = pow(finalColor, vec3(0.95));
    
    FragColor = vec4(finalColor, 1.0);
}
