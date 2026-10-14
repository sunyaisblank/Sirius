#version 330 core
in vec2 TexCoord;
out vec4 FragColor;

uniform sampler2D screenTexture;

// The texture contains display-linear sRGB after the host-selected exposure,
// tone map, colour grade, and optional film finish. This shader owns only the
// final IEC 61966-2-1 transfer encode for the live window. OpenGL framebuffer
// sRGB conversion is explicitly disabled by the caller, so this is the sole
// live-view transfer boundary.
const float kSrgbLinearBreakpoint = 0.0031308;
const float kSrgbLinearSlope = 12.92;
const float kSrgbPowerScale = 1.055;
const float kSrgbPowerOffset = 0.055;
const float kSrgbPowerExponent = 1.0 / 2.4;

float EncodeSrgbChannel(float linear) {
    linear = clamp(linear, 0.0, 1.0);
    if (linear <= kSrgbLinearBreakpoint) {
        return kSrgbLinearSlope * linear;
    }
    return kSrgbPowerScale * pow(linear, kSrgbPowerExponent) - kSrgbPowerOffset;
}

void main() {
    vec3 displayLinear = texture(screenTexture, TexCoord).rgb;
    vec3 encoded = vec3(
        EncodeSrgbChannel(displayLinear.r),
        EncodeSrgbChannel(displayLinear.g),
        EncodeSrgbChannel(displayLinear.b));
    FragColor = vec4(encoded, 1.0);
}
