// RDOP006A.h - Unified OptiX C API Declarations
// Component ID: RDOP006A
// Single header for all extern "C" OptiX host-side entry points.
// Included by both RDRT001A.cpp (Renderer) and RDOX001A.h (OptiXAccelerator).

#pragma once

#include "RDOP003A.h"

extern "C" {
    typedef void* SiriusOptixHandle;

    // Lifecycle
    SiriusOptixHandle sirius_optix_create();
    void sirius_optix_destroy(SiriusOptixHandle handle);
    bool sirius_optix_initialize(SiriusOptixHandle handle, int width, int height);
    bool sirius_optix_create_pipeline(SiriusOptixHandle handle, const char* ptxPath);
    void sirius_optix_cleanup(SiriusOptixHandle handle);
    bool sirius_optix_is_initialized(SiriusOptixHandle handle);

    // Rendering
    void sirius_optix_launch(SiriusOptixHandle handle, const Sirius::LaunchParams* params);
    void sirius_optix_update_display(SiriusOptixHandle handle);
    void sirius_optix_resize(SiriusOptixHandle handle, int width, int height);
    void sirius_optix_set_metric_type(SiriusOptixHandle handle, int type);

    // Framebuffer
    float* sirius_optix_get_frame_buffer(SiriusOptixHandle handle);
    void sirius_optix_register_gl_texture(SiriusOptixHandle handle, unsigned int glTexture);

    // Background texture
    bool sirius_optix_upload_background(SiriusOptixHandle handle, const unsigned char* data, int width, int height);
    unsigned long long sirius_optix_get_background_texture(SiriusOptixHandle handle);

    // Denoiser
    bool sirius_optix_init_denoiser(SiriusOptixHandle handle);
    void sirius_optix_denoise(SiriusOptixHandle handle, float blendFactor);
    void sirius_optix_set_denoiser_enabled(SiriusOptixHandle handle, bool enabled);
}
