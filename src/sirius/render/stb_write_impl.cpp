// Single translation unit for the stb_image_write encoder implementation.
// Keeping vendored implementation code separate lets the first-party PNG
// writer retain strict warnings on every supported compiler.

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
