// End-to-end CPU render-session probe (session-level gate).
//
// Drives RenderSession CPU-only at 64x64, 4 spp, Kerr a=0.9, writing a PNG and
// an EXR to a temp directory. Asserts: the session reaches Complete, both files
// exist, the PNG decodes (stb) and the EXR loads (tinyexr), and each decoded
// image is finite and non-constant. Full pixel-identity against the reference
// tapes happens at the APP stage with the real CLI; this probe pins that the
// ported session, tracer, and writers wire together and produce a real image.

#include "sirius/render/session/render_session.h"

#include <gtest/gtest.h>

#include <stb_image.h>
#include <tinyexr.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <string>

namespace {

using sirius::render::RenderSession;
using sirius::render::SessionConfig;
using sirius::render::SessionState;

SessionConfig ProbeConfig(const std::string& output_path) {
    SessionConfig cfg;
    cfg.width = 64;
    cfg.height = 64;
    cfg.samplesPerPixel = 4;
    cfg.tileSize = 64;
    cfg.enableParallelRendering = false;  // Deterministic single-thread for the probe.
    cfg.metricId = sirius::core::MetricId::Kerr;
    cfg.blackHoleMass = 1.0;
    cfg.blackHoleSpin = 0.9;
    cfg.outputPath = output_path;
    return cfg;
}

// Run one CPU render to output_path; returns the terminal session state.
SessionState RenderTo(const std::string& output_path) {
    RenderSession session;
    if (!session.Configure(ProbeConfig(output_path))) {
        return SessionState::Failed;
    }
    return session.Execute();
}

}  // namespace

TEST(RenderSessionProbe, CpuKerrRenderProducesValidPngAndExr) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "sirius_render_probe";
    fs::create_directories(dir);

    const std::string pngPath = (dir / "probe_kerr.png").string();
    const std::string exrPath = (dir / "probe_kerr.exr").string();
    fs::remove(pngPath);
    fs::remove(exrPath);

    // --- PNG render ---------------------------------------------------------
    ASSERT_EQ(RenderTo(pngPath), SessionState::Complete) << "PNG render did not complete";
    ASSERT_TRUE(fs::exists(pngPath)) << "PNG file was not written";

    int pw = 0, ph = 0, pc = 0;
    unsigned char* png = stbi_load(pngPath.c_str(), &pw, &ph, &pc, 3);
    ASSERT_NE(png, nullptr) << "PNG failed to decode: " << stbi_failure_reason();
    EXPECT_EQ(pw, 64);
    EXPECT_EQ(ph, 64);

    // PNG (8-bit) is always finite; require it to be non-constant.
    bool png_varies = false;
    for (int i = 3; i < pw * ph * 3 && !png_varies; ++i) {
        if (png[i] != png[0]) png_varies = true;
    }
    EXPECT_TRUE(png_varies) << "PNG image is constant (nothing rendered)";
    stbi_image_free(png);

    // --- EXR render ---------------------------------------------------------
    ASSERT_EQ(RenderTo(exrPath), SessionState::Complete) << "EXR render did not complete";
    ASSERT_TRUE(fs::exists(exrPath)) << "EXR file was not written";

    float* exr = nullptr;
    int ew = 0, eh = 0;
    const char* err = nullptr;
    int ret = LoadEXR(&exr, &ew, &eh, exrPath.c_str(), &err);
    ASSERT_EQ(ret, TINYEXR_SUCCESS) << (err ? err : "unknown tinyexr error");
    EXPECT_EQ(ew, 64);
    EXPECT_EQ(eh, 64);

    // EXR carries linear HDR floats; require every sample finite and the frame
    // non-constant.
    bool exr_finite = true;
    bool exr_varies = false;
    const float first = exr[0];
    for (int i = 0; i < ew * eh * 4; ++i) {
        if (!std::isfinite(exr[i])) exr_finite = false;
        if ((i % 4) != 3 && exr[i] != first) exr_varies = true;  // Ignore the alpha channel.
    }
    EXPECT_TRUE(exr_finite) << "EXR contains non-finite samples";
    EXPECT_TRUE(exr_varies) << "EXR image is constant (nothing rendered)";
    std::free(exr);

    fs::remove(pngPath);
    fs::remove(exrPath);
}
