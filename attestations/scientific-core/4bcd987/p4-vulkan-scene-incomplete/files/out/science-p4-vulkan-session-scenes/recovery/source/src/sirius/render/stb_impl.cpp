// Single translation unit for the stb_image decoder. Ported from RDST001A.cpp.
//
// The session's starfield loader (and the PNG round-trip tests) need the
// decoder; defining STB_IMAGE_IMPLEMENTATION exactly once here keeps its
// symbols out of every other TU. The encoder side lives in stb_write_impl.cpp.

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
