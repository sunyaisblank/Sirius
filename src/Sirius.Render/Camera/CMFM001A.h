// CMFM001A.h - Film Format Definitions
// Component ID: CMFM001A
// Physical film formats for camera resolution calculation
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_CMFM001A_H
#define SIRIUS_RENDER_CMFM001A_H

#include <cmath>
#include <string>

namespace sirius::render {

//==============================================================================
// Resolution
// Resolution structure for computed dimensions
//==============================================================================

struct Resolution {
    int width = 0;
    int height = 0;
    double megapixels = 0.0;

    Resolution() = default;
    Resolution(int w, int h)
        : width(w), height(h)
        , megapixels(static_cast<double>(w) * h / 1e6) {}
};

//==============================================================================
// FilmFormat
// Physical film gate dimensions
//==============================================================================

struct FilmFormat {
    std::string name;
    double gateWidth = 36.0;    // Gate width in mm
    double gateHeight = 24.0;   // Gate height in mm
    double aspectRatio = 1.5;   // Width / Height

    FilmFormat() {
        aspectRatio = gateWidth / gateHeight;
    }

    FilmFormat(const std::string& n, double w, double h)
        : name(n), gateWidth(w), gateHeight(h)
        , aspectRatio(w / h) {}

    //--------------------------------------------------------------------------
    // Standard Film Formats
    //--------------------------------------------------------------------------

    /// IMAX 15-perf 70mm (horizontal)
    static FilmFormat IMAX70mm() {
        return FilmFormat("IMAX 70mm 15/70", 69.6, 48.5);
    }

    /// IMAX GT 15-perf 70mm (dome format)
    static FilmFormat IMAXGT() {
        return FilmFormat("IMAX GT", 69.6, 48.5);
    }

    /// VistaVision 8-perf 35mm (horizontal)
    static FilmFormat VistaVision() {
        return FilmFormat("VistaVision", 37.72, 25.17);
    }

    /// Standard 35mm Full Frame
    static FilmFormat FullFrame35mm() {
        return FilmFormat("35mm Full Frame", 36.0, 24.0);
    }

    /// Super 35mm
    static FilmFormat Super35mm() {
        return FilmFormat("Super 35mm", 24.89, 18.66);
    }

    /// ARRI Alexa 65 sensor
    static FilmFormat Alexa65() {
        return FilmFormat("ARRI Alexa 65", 54.12, 25.58);
    }

    /// RED Monstro 8K VV
    static FilmFormat REDMonstro() {
        return FilmFormat("RED Monstro 8K VV", 40.96, 21.60);
    }

    /// Hasselblad H6D-100c (medium format)
    static FilmFormat MediumFormat() {
        return FilmFormat("Medium Format 645", 53.4, 40.0);
    }

    /// Large Format 4x5
    static FilmFormat LargeFormat4x5() {
        return FilmFormat("4x5 Large Format", 121.0, 97.0);
    }

    /// 8x10 Large Format
    static FilmFormat LargeFormat8x10() {
        return FilmFormat("8x10 Large Format", 254.0, 203.0);
    }
};

//==============================================================================
// Camera Presets
// Pre-configured camera setups for common use cases
//==============================================================================

namespace presets {

/// Standard IMAX 70mm camera with 50mm lens
struct IMAXCamera {
    FilmFormat format;
    double focalLength;         // Focal length in mm
    std::string name;

    double horizontalFOV() const {
        return 2.0 * std::atan(format.gateWidth / (2.0 * focalLength)) * 180.0 / M_PI;
    }

    double verticalFOV() const {
        return 2.0 * std::atan(format.gateHeight / (2.0 * focalLength)) * 180.0 / M_PI;
    }

    Resolution calculateResolution(double pixelsPerMM) const {
        int w = static_cast<int>(format.gateWidth * pixelsPerMM);
        int h = static_cast<int>(format.gateHeight * pixelsPerMM);
        return Resolution(w, h);
    }
};

/// IMAX 70mm with standard 50mm lens
inline IMAXCamera IMAX_70mm_Standard() {
    return IMAXCamera{FilmFormat::IMAX70mm(), 50.0, "IMAX 70mm Standard"};
}

/// IMAX 70mm with wide 35mm lens
inline IMAXCamera IMAX_70mm_Wide() {
    return IMAXCamera{FilmFormat::IMAX70mm(), 35.0, "IMAX 70mm Wide"};
}

/// VistaVision with 50mm lens
inline IMAXCamera VistaVision_Standard() {
    return IMAXCamera{FilmFormat::VistaVision(), 50.0, "VistaVision Standard"};
}

/// Alexa 65 with 50mm lens
inline IMAXCamera Alexa65_Standard() {
    return IMAXCamera{FilmFormat::Alexa65(), 50.0, "ARRI Alexa 65 Standard"};
}

} // namespace presets

} // namespace sirius::render

#endif // SIRIUS_RENDER_CMFM001A_H
