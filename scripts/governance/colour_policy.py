"""Single-authority wiring for spectral radiance and display encoding."""

from __future__ import annotations

import re
from pathlib import Path

from .common import (
    CPP_NON_CODE,
    SOURCE_ROOT,
    relative,
)


SRGB_TRANSFER_AUTHORITY = SOURCE_ROOT / "core" / "srgb_transfer.h"
SRGB_TRANSFER_VIEWER_MIRROR = (
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD003A.frag"
)
VIEWER_COMMAND = SOURCE_ROOT / "app" / "cli" / "view_command.cpp"
LIVE_VIEWER_SHADERS = {
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD003A.frag",
    SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD003A.vert",
}
SRGB_TRANSFER_CONSUMERS = {
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h": "colour::EncodeSrgbChannel",
    SOURCE_ROOT / "render" / "image_buffer.h": "core::colour::TryEncodeSrgb8",
    SOURCE_ROOT / "render" / "png_writer.h": "core::colour::TryEncodeSrgb8",
    SOURCE_ROOT / "render" / "exr_writer.h": "core::colour::TryEncodeSrgb8",
}
SRGB_LINEAR_BREAKPOINT = re.compile(
    r"(?<![A-Za-z0-9_.])0\.0031308(?:[fFlL])?(?![A-Za-z0-9_.])"
)
APPROXIMATE_GAMMA_22 = re.compile(r"\bpow\s*\([^;]*1\.0\s*/\s*2\.2")
FIXED_REINHARD_ASSIGNMENT = re.compile(
    r"\b(?P<value>[A-Za-z_]\w*)\s*=\s*(?P=value)\s*/\s*"
    r"\(\s*(?P=value)\s*\+\s*vec3\s*\(\s*1\.0\s*\)\s*\)"
)
XYZ_SRGB_AUTHORITY = SOURCE_ROOT / "core" / "xyz_srgb.h"
XYZ_SRGB_HOST_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
)
XYZ_SRGB_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_disk.slang"
XYZ_SRGB_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
XYZ_SRGB_COEFFICIENT_PAIRS = (
    ("12831", "3959"),
    ("329", "214"),
    ("1974", "3959"),
    ("851781", "878810"),
    ("1648619", "878810"),
    ("36519", "878810"),
    ("705", "12673"),
    ("2585", "12673"),
    ("705", "667"),
)
XYZ_SRGB_LEGACY_COEFFICIENT = re.compile(
    r"(?<![A-Za-z0-9_.])(?:3\.2404542|1\.5371385|0\.4985314|0\.9692660|"
    r"1\.8760108|0\.0415560|0\.0556434|0\.2040259|1\.0572252)"
    r"(?:[fFlL])?(?![A-Za-z0-9_.])"
)
CIE1931_OBSERVER_AUTHORITY = SOURCE_ROOT / "core" / "cie1931_observer.h"
CIE1931_OBSERVER_HOST_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
    SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h",
)
CIE1931_OBSERVER_DEVICE_AUTHORITY = SOURCE_ROOT / "kernels" / "gr_disk.slang"
CIE1931_OBSERVER_PARITY_PROBE = SOURCE_ROOT / "kernels" / "parity_probe.slang"
CIE1931_GAUSSIAN_LOBES = (
    ("442.0", "0.0624", "0.0374", "0.362"),
    ("599.8", "0.0264", "0.0323", "1.056"),
    ("501.1", "0.0490", "0.0382", "0.065"),
    ("568.8", "0.0213", "0.0247", "0.821"),
    ("530.9", "0.0613", "0.0322", "0.286"),
    ("437.0", "0.0845", "0.0278", "1.217"),
    ("459.0", "0.0385", "0.0725", "0.681"),
)
CIE1931_LEGACY_FACADE = re.compile(r"\bCie[XYZ]\s*\(")
CIE1931_STALE_TABLE = re.compile(
    r"(?<![A-Za-z0-9_.])0\.014310(?:[fFlL])?.{0,160}"
    r"0\.043510(?:[fFlL])?.{0,160}0\.134380(?:[fFlL])?",
    re.DOTALL,
)
BLACKBODY_LAWS_AUTHORITY = SOURCE_ROOT / "core" / "spectral" / "blackbody_laws.h"
BLACKBODY_LAWS_CONSUMERS = (
    SOURCE_ROOT / "core" / "spectral" / "blackbody.h",
    SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h",
)
PHYSICAL_CONSTANTS_AUTHORITY = SOURCE_ROOT / "core" / "constants.h"
PLANCK_EXPM1 = re.compile(r"\bstd::expm1\s*\(\s*exponent\s*\)")
PLANCK_CONSTANT_IDENTIFIER = re.compile(r"\bk(?:Planck(?:C[12])?|Boltzmann)\b")
WIEN_CONSTANT_LITERAL = re.compile(
    r"(?<![A-Za-z0-9_.])2\.897771955e-3(?:[fFlL])?(?![A-Za-z0-9_.])",
    re.IGNORECASE,
)
ACES_FIT_AUTHORITY = SOURCE_ROOT / "core" / "postprocess.h"
ACES_FIT_CONSUMERS = {
    SOURCE_ROOT / "app" / "config" / "config_schema.h": "kDefaultTonemapperName",
    SOURCE_ROOT / "app" / "config" / "config_loader.cpp": "ParseTonemapType",
    SOURCE_ROOT / "app" / "config" / "session_config_adapter.cpp": "ParseTonemapType",
    SOURCE_ROOT / "render" / "session" / "render_session.h": "TonemapType::AcesFit",
    SOURCE_ROOT / "render" / "session" / "render_session.cpp": "TonemapType::AcesFit",
}
ACES_FIT_CLI = SOURCE_ROOT / "app" / "cli" / "render_command.cpp"
SPECTRAL_RADIANCE_FACADE = SOURCE_ROOT / "core" / "spectral" / "spectral_radiance.h"
ACES_FIT_COEFFICIENTS = ("2.51", "0.03", "2.43", "0.59", "0.14")
LEGACY_ACES_TONEMAP_IDENTIFIER = re.compile(
    r"\bTonemapType::Aces\b|\btonemap::Aces\s*\("
)
FALSE_SPECTRAL_ACES_API = re.compile(r"\bstruct\s+Aces\b|\bToAces\s*\(")
FALSE_ABSOLUTE_SPECTRAL_SRGB_API = re.compile(r"\bToSrgb\s*\(")
ABSOLUTE_SPECTRAL_DISPLAY_ENCODING = re.compile(
    r"\b(?:XyzD65ToLinearSrgb|EncodeSrgbChannel|EncodeClippedSrgbChannel|TryEncodeSrgb8)\s*\("
)
DETACHED_XYZ_AP0_MATRIX = re.compile(
    r"1\.0498110175.{0,240}-0\.4959030231.{0,240}"
    r"1\.3733130458.{0,240}0\.9912520182",
    re.DOTALL,
)
BARE_ACES_CONFIG_LITERAL = re.compile(r'"ACES"')


def srgb_transfer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(SRGB_TRANSFER_AUTHORITY)
    if authority is None:
        return ["the host sRGB transfer authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "EncodeSrgbChannel",
        "EncodeClippedSrgbChannel",
        "TryEncodeSrgb8",
    ):
        if marker not in authority_code:
            errors.append(f"the host sRGB transfer authority omits {marker}")
    if len(SRGB_LINEAR_BREAKPOINT.findall(authority_code)) != 1:
        errors.append("the host sRGB transfer authority must own one IEC linear breakpoint")

    viewer = documents.get(SRGB_TRANSFER_VIEWER_MIRROR)
    if viewer is None:
        errors.append("the live viewer sRGB transfer mirror is missing")
    else:
        viewer_code = CPP_NON_CODE.sub(" ", viewer)
        for marker in (
            "const float kSrgbLinearBreakpoint = 0.0031308;",
            "const float kSrgbLinearSlope = 12.92;",
            "const float kSrgbPowerScale = 1.055;",
            "const float kSrgbPowerOffset = 0.055;",
            "const float kSrgbPowerExponent = 1.0 / 2.4;",
            "EncodeSrgbChannel",
            "linear = clamp(linear, 0.0, 1.0);",
            "if (linear <= kSrgbLinearBreakpoint)",
            "return kSrgbLinearSlope * linear;",
            "kSrgbPowerScale * pow(linear, kSrgbPowerExponent) - kSrgbPowerOffset",
            "vec3 displayLinear = texture(screenTexture, TexCoord).rgb;",
            "FragColor = vec4(encoded, 1.0);",
        ):
            if marker not in viewer_code:
                errors.append(f"the live viewer transfer mirror omits {marker}")
        if len(SRGB_LINEAR_BREAKPOINT.findall(viewer_code)) != 1:
            errors.append("the live viewer mirror must own one IEC linear breakpoint")
        if viewer_code.count("EncodeSrgbChannel(displayLinear.") != 3:
            errors.append("the live viewer must transfer-encode exactly three RGB channels")
        if viewer_code.count("texture(") != 1:
            errors.append("the live viewer transfer shader must sample exactly one texture value")
        if viewer_code.count("pow(") != 1:
            errors.append("the live viewer transfer shader must own one IEC power branch")
        if APPROXIMATE_GAMMA_22.search(viewer_code):
            errors.append("the live viewer restored an approximate gamma-2.2 encode")
        if FIXED_REINHARD_ASSIGNMENT.search(viewer_code):
            errors.append("the live viewer restored an unconfigured secondary tone map")

    view_command = documents.get(VIEWER_COMMAND)
    if view_command is None:
        errors.append("the live viewer command is missing")
    else:
        command_code = CPP_NON_CODE.sub(" ", view_command)
        if command_code.count("glDisable(GL_FRAMEBUFFER_SRGB)") != 1:
            errors.append("the live viewer must explicitly disable hardware sRGB re-encoding")
        if "glEnable(GL_FRAMEBUFFER_SRGB)" in command_code:
            errors.append("the live viewer enables a second hardware sRGB encode")
        if command_code.count("GL_RGBA32F") != 1:
            errors.append("the live viewer must preserve display-linear values in one fp32 texture")

    required_include = '#include "sirius/core/srgb_transfer.h"'
    for path, marker in SRGB_TRANSFER_CONSUMERS.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"sRGB transfer consumer is missing: {relative(path)}")
            continue
        if required_include not in document or marker not in CPP_NON_CODE.sub(" ", document):
            errors.append(
                f"{relative(path)} does not delegate to the host sRGB transfer authority"
            )

    spectral_radiance = documents.get(SPECTRAL_RADIANCE_FACADE)
    if spectral_radiance is None:
        errors.append("the binned physical spectral-radiance facade is missing")
    else:
        spectral_code = CPP_NON_CODE.sub(" ", spectral_radiance)
        if FALSE_ABSOLUTE_SPECTRAL_SRGB_API.search(spectral_code):
            errors.append(
                "the binned physical spectral-radiance facade exposes an sRGB display API"
            )
        if ABSOLUTE_SPECTRAL_DISPLAY_ENCODING.search(spectral_code):
            errors.append(
                "the binned physical spectral-radiance facade directly performs display encoding"
            )

    for path, document in documents.items():
        if path in {SRGB_TRANSFER_AUTHORITY, SRGB_TRANSFER_VIEWER_MIRROR}:
            continue
        if SRGB_LINEAR_BREAKPOINT.search(CPP_NON_CODE.sub(" ", document)):
            errors.append(f"{relative(path)} reimplements the IEC sRGB transfer breakpoint")
    return errors


def verify_srgb_transfer_authority_policy() -> None:
    valid = {
        SRGB_TRANSFER_AUTHORITY: (
            "EncodeSrgbChannel EncodeClippedSrgbChannel TryEncodeSrgb8 0.0031308"
        ),
        **{
            path: '#include "sirius/core/srgb_transfer.h"\n' + marker
            for path, marker in SRGB_TRANSFER_CONSUMERS.items()
        },
        SPECTRAL_RADIANCE_FACADE: "struct SpectralRadiance { Xyz ToXyz() const; };",
        SRGB_TRANSFER_VIEWER_MIRROR: (
            "const float kSrgbLinearBreakpoint = 0.0031308;\n"
            "const float kSrgbLinearSlope = 12.92;\n"
            "const float kSrgbPowerScale = 1.055;\n"
            "const float kSrgbPowerOffset = 0.055;\n"
            "const float kSrgbPowerExponent = 1.0 / 2.4;\n"
            "float EncodeSrgbChannel(float linear) {\n"
            " linear = clamp(linear, 0.0, 1.0);\n"
            " if (linear <= kSrgbLinearBreakpoint)\n"
            "  return kSrgbLinearSlope * linear;\n"
            " return kSrgbPowerScale * pow(linear, kSrgbPowerExponent) - "
            "kSrgbPowerOffset;\n}\n"
            "void main() { vec3 displayLinear = texture(screenTexture, TexCoord).rgb;\n"
            " vec3 encoded = vec3(EncodeSrgbChannel(displayLinear.r), "
            "EncodeSrgbChannel(displayLinear.g), EncodeSrgbChannel(displayLinear.b));\n"
            " FragColor = vec4(encoded, 1.0); }\n"
        ),
        VIEWER_COMMAND: "glDisable(GL_FRAMEBUFFER_SRGB); GL_RGBA32F",
    }
    if srgb_transfer_authority_errors(valid):
        raise RuntimeError("sRGB transfer policy rejected the single host authority")
    duplicated = dict(valid)
    consumer = next(iter(SRGB_TRANSFER_CONSUMERS))
    duplicated[consumer] += "\n0.0031308"
    if not srgb_transfer_authority_errors(duplicated):
        raise RuntimeError("sRGB transfer policy accepted a production reimplementation")
    detached = dict(valid)
    detached[consumer] = "independent_transfer();"
    if not srgb_transfer_authority_errors(detached):
        raise RuntimeError("sRGB transfer policy accepted a detached production consumer")

    false_spectral_api = dict(valid)
    false_spectral_api[SPECTRAL_RADIANCE_FACADE] += "\nRgb ToSrgb() const;"
    if not srgb_transfer_authority_errors(false_spectral_api):
        raise RuntimeError("sRGB policy accepted an absolute-radiance display facade")

    direct_spectral_encoding = dict(valid)
    direct_spectral_encoding[SPECTRAL_RADIANCE_FACADE] += (
        "\ncolour::XyzD65ToLinearSrgb(x, y, z); colour::EncodeSrgbChannel(r);"
    )
    if not srgb_transfer_authority_errors(direct_spectral_encoding):
        raise RuntimeError("sRGB policy accepted direct encoding of absolute spectral radiance")

    double_tonemap = dict(valid)
    double_tonemap[SRGB_TRANSFER_VIEWER_MIRROR] += (
        "\ncolor = color / (color + vec3(1.0));\n"
    )
    if not srgb_transfer_authority_errors(double_tonemap):
        raise RuntimeError("sRGB policy accepted a secondary viewer tone map")

    double_encode = dict(valid)
    double_encode[VIEWER_COMMAND] += "\nglEnable(GL_FRAMEBUFFER_SRGB);\n"
    if not srgb_transfer_authority_errors(double_encode):
        raise RuntimeError("sRGB policy accepted double viewer transfer encoding")

    if set(LIVE_VIEWER_SHADERS) != {
            SRGB_TRANSFER_VIEWER_MIRROR,
            SOURCE_ROOT / "app" / "viewer" / "shaders" / "RDSD003A.vert"}:
        raise RuntimeError("live viewer shader inventory policy is internally inconsistent")


def xyz_srgb_ratio(numerator: str, denominator: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?<![A-Za-z0-9_.]){numerator}(?:\.0)?(?:[fFlL])?\s*/\s*"
        rf"{denominator}(?:\.0)?(?:[fFlL])?(?![A-Za-z0-9_.])"
    )


def xyz_srgb_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(XYZ_SRGB_AUTHORITY)
    if authority is None:
        return ["the host XYZ-D65 to linear-sRGB authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    if "XyzD65ToLinearSrgb" not in authority_code:
        errors.append("the host XYZ-D65 to linear-sRGB authority omits its transform")
    for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
        if len(xyz_srgb_ratio(numerator, denominator).findall(authority_code)) != 1:
            errors.append(
                "the host XYZ-D65 to linear-sRGB authority must own exact coefficient "
                f"{numerator}/{denominator} once"
            )

    required_include = '#include "sirius/core/xyz_srgb.h"'
    for path in XYZ_SRGB_HOST_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"XYZ-to-sRGB host consumer is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        if required_include not in document or "colour::XyzD65ToLinearSrgb" not in code:
            errors.append(
                f"{relative(path)} does not delegate to the host XYZ-to-sRGB authority"
            )

    device = documents.get(XYZ_SRGB_DEVICE_AUTHORITY)
    if device is None:
        errors.append("the Slang XYZ-D65 to linear-sRGB mirror is missing")
    else:
        device_code = CPP_NON_CODE.sub(" ", device)
        if device_code.count("XyzD65ToLinearSrgb") < 2:
            errors.append("the Slang blackbody path does not consume its XYZ-to-sRGB mirror")
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
            if len(xyz_srgb_ratio(numerator, denominator).findall(device_code)) != 1:
                errors.append(
                    "the Slang XYZ-D65 to linear-sRGB mirror must own exact coefficient "
                    f"{numerator}/{denominator} once"
                )

    parity = documents.get(XYZ_SRGB_PARITY_PROBE)
    if parity is None or "OP_XYZ_D65_TO_LINEAR_SRGB" not in parity or (
        "XyzD65ToLinearSrgb" not in CPP_NON_CODE.sub(" ", parity)
    ):
        errors.append("the direct XYZ-D65 to linear-sRGB device parity route is missing")

    allowed_ratio_paths = {XYZ_SRGB_AUTHORITY, XYZ_SRGB_DEVICE_AUTHORITY}
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if XYZ_SRGB_LEGACY_COEFFICIENT.search(code):
            errors.append(f"{relative(path)} retains a rounded legacy XYZ-to-sRGB coefficient")
        if path in allowed_ratio_paths:
            continue
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS:
            if xyz_srgb_ratio(numerator, denominator).search(code):
                errors.append(
                    f"{relative(path)} reimplements exact XYZ-to-sRGB coefficient "
                    f"{numerator}/{denominator}"
                )
    return errors


def verify_xyz_srgb_authority_policy() -> None:
    host_coefficients = " ".join(
        f"{numerator}.0L / {denominator}.0L"
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS
    )
    device_coefficients = " ".join(
        f"{numerator}.0f / {denominator}.0f"
        for numerator, denominator in XYZ_SRGB_COEFFICIENT_PAIRS
    )
    valid = {
        XYZ_SRGB_AUTHORITY: "XyzD65ToLinearSrgb " + host_coefficients,
        **{
            path: '#include "sirius/core/xyz_srgb.h"\ncolour::XyzD65ToLinearSrgb'
            for path in XYZ_SRGB_HOST_CONSUMERS
        },
        XYZ_SRGB_DEVICE_AUTHORITY: (
            "public float3 XyzD65ToLinearSrgb "
            + device_coefficients
            + " BlackbodyColor XyzD65ToLinearSrgb"
        ),
        XYZ_SRGB_PARITY_PROBE: "OP_XYZ_D65_TO_LINEAR_SRGB XyzD65ToLinearSrgb",
    }
    if xyz_srgb_authority_errors(valid):
        raise RuntimeError("XYZ-to-sRGB policy rejected the governed host/device authority")

    duplicated = dict(valid)
    consumer = XYZ_SRGB_HOST_CONSUMERS[0]
    duplicated[consumer] += "\n12831.0 / 3959.0"
    if not xyz_srgb_authority_errors(duplicated):
        raise RuntimeError("XYZ-to-sRGB policy accepted a production coefficient copy")

    legacy = dict(valid)
    legacy[consumer] += "\n3.2404542f"
    if not xyz_srgb_authority_errors(legacy):
        raise RuntimeError("XYZ-to-sRGB policy accepted a rounded legacy matrix")

    detached = dict(valid)
    detached[consumer] = "independent_xyz_to_rgb();"
    if not xyz_srgb_authority_errors(detached):
        raise RuntimeError("XYZ-to-sRGB policy accepted a detached host consumer")


def cie1931_lobe(values: tuple[str, str, str, str]) -> re.Pattern[str]:
    def literal(value: str) -> str:
        whole, dot, fraction = value.partition(".")
        decimal = rf"{whole}\.{fraction}" if dot else whole
        return rf"(?<![A-Za-z0-9_.]){decimal}(?:[fFlL])?(?![A-Za-z0-9_.])"

    return re.compile(r"\s*,\s*".join(literal(value) for value in values))


def cie1931_observer_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(CIE1931_OBSERVER_AUTHORITY)
    if authority is None:
        return ["the host CIE 1931 observer-fit authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in ("CieXyzMatching", "Cie1931TwoDegreeFit", "std::isfinite"):
        if marker not in authority_code:
            errors.append(f"the host CIE 1931 observer-fit authority omits {marker}")
    for lobe in CIE1931_GAUSSIAN_LOBES:
        if len(cie1931_lobe(lobe).findall(authority_code)) != 1:
            errors.append(
                "the host CIE 1931 observer-fit authority must own Gaussian lobe "
                + "/".join(lobe)
                + " once"
            )

    required_include = '#include "sirius/core/cie1931_observer.h"'
    for path in CIE1931_OBSERVER_HOST_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"CIE 1931 observer-fit host consumer is missing: {relative(path)}")
            continue
        code = CPP_NON_CODE.sub(" ", document)
        if required_include not in document or "colour::Cie1931TwoDegreeFit" not in code:
            errors.append(
                f"{relative(path)} does not delegate to the host CIE 1931 observer-fit authority"
            )

    device = documents.get(CIE1931_OBSERVER_DEVICE_AUTHORITY)
    if device is None:
        errors.append("the Slang CIE 1931 observer-fit mirror is missing")
    else:
        device_code = CPP_NON_CODE.sub(" ", device)
        if device_code.count("Cie1931TwoDegreeFit") < 2:
            errors.append("the Slang blackbody path does not consume its CIE observer-fit mirror")
        for lobe in CIE1931_GAUSSIAN_LOBES:
            if len(cie1931_lobe(lobe).findall(device_code)) != 1:
                errors.append(
                    "the Slang CIE 1931 observer-fit mirror must own Gaussian lobe "
                    + "/".join(lobe)
                    + " once"
                )

    parity = documents.get(CIE1931_OBSERVER_PARITY_PROBE)
    if parity is None or "OP_CIE_1931_TWO_DEGREE_FIT" not in parity or (
        "Cie1931TwoDegreeFit" not in CPP_NON_CODE.sub(" ", parity)
    ):
        errors.append("the direct CIE 1931 observer-fit device parity route is missing")

    allowed_lobe_paths = {
        CIE1931_OBSERVER_AUTHORITY,
        CIE1931_OBSERVER_DEVICE_AUTHORITY,
    }
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if CIE1931_LEGACY_FACADE.search(code):
            errors.append(f"{relative(path)} retains a detached CIE colour-matching facade")
        if CIE1931_STALE_TABLE.search(code) or re.search(r"\bk[XYZ]Bar\b", code):
            errors.append(f"{relative(path)} retains the stale wavelength-mislabeled CIE table")
        if path in allowed_lobe_paths:
            continue
        for lobe in CIE1931_GAUSSIAN_LOBES:
            if cie1931_lobe(lobe).search(code):
                errors.append(
                    f"{relative(path)} reimplements CIE observer-fit Gaussian lobe "
                    + "/".join(lobe)
                )
    return errors


def verify_cie1931_observer_authority_policy() -> None:
    coefficients = " ".join(", ".join(lobe) for lobe in CIE1931_GAUSSIAN_LOBES)
    valid = {
        CIE1931_OBSERVER_AUTHORITY: (
            "CieXyzMatching Cie1931TwoDegreeFit std::isfinite " + coefficients
        ),
        **{
            path: '#include "sirius/core/cie1931_observer.h"\ncolour::Cie1931TwoDegreeFit'
            for path in CIE1931_OBSERVER_HOST_CONSUMERS
        },
        CIE1931_OBSERVER_DEVICE_AUTHORITY: (
            "public float3 Cie1931TwoDegreeFit "
            + coefficients
            + " BlackbodyColor Cie1931TwoDegreeFit"
        ),
        CIE1931_OBSERVER_PARITY_PROBE: (
            "OP_CIE_1931_TWO_DEGREE_FIT Cie1931TwoDegreeFit"
        ),
    }
    if cie1931_observer_authority_errors(valid):
        raise RuntimeError("CIE observer policy rejected the governed host/device authority")

    consumer = CIE1931_OBSERVER_HOST_CONSUMERS[0]
    duplicated = dict(valid)
    duplicated[consumer] += "\n442.0, 0.0624, 0.0374, 0.362"
    if not cie1931_observer_authority_errors(duplicated):
        raise RuntimeError("CIE observer policy accepted a production coefficient copy")

    detached = dict(valid)
    detached[consumer] = "independent_cie_observer();"
    if not cie1931_observer_authority_errors(detached):
        raise RuntimeError("CIE observer policy accepted a detached host consumer")

    stale_table = dict(valid)
    stale_table[consumer] += "\n0.014310, 0.043510, 0.134380"
    if not cie1931_observer_authority_errors(stale_table):
        raise RuntimeError("CIE observer policy accepted the wavelength-mislabeled table")

    missing_parity = dict(valid)
    missing_parity[CIE1931_OBSERVER_PARITY_PROBE] = "unrelated_probe();"
    if not cie1931_observer_authority_errors(missing_parity):
        raise RuntimeError("CIE observer policy accepted a missing direct device parity route")


def aces_fit_sequence() -> re.Pattern[str]:
    literals = [
        rf"(?<![A-Za-z0-9_.]){re.escape(value)}(?:[fFlL])?(?![A-Za-z0-9_.])"
        for value in ACES_FIT_COEFFICIENTS
    ]
    return re.compile(r".{0,180}".join(literals), re.DOTALL)


def aces_contract_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(ACES_FIT_AUTHORITY)
    if authority is None:
        return ["the host ACES-fit tone-map authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "AcesFit",
        "ParseTonemapType",
        "kDefaultTonemapperName",
        "kTonemapNames",
        "SupportedTonemapperNames",
        "std::isfinite",
    ):
        if marker not in authority_code:
            errors.append(f"the host ACES-fit tone-map authority omits {marker}")
    if len(aces_fit_sequence().findall(authority_code)) != 1:
        errors.append("the host ACES-fit authority must own the five fit coefficients once")
    if "knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve" not in authority:
        errors.append("the ACES-fit authority does not identify the original fit source")

    for path, marker in ACES_FIT_CONSUMERS.items():
        document = documents.get(path)
        if document is None:
            errors.append(f"ACES-fit consumer is missing: {relative(path)}")
            continue
        if marker not in CPP_NON_CODE.sub(" ", document):
            errors.append(f"{relative(path)} does not consume the explicit ACES-fit authority")

    cli = documents.get(ACES_FIT_CLI)
    if cli is None:
        errors.append("the render CLI tone-map help is missing")
    elif "Tonemapper: ACESFit" not in cli or "not an ACES Output Transform" not in cli:
        errors.append("the render CLI does not distinguish ACESFit from an ACES Output Transform")

    bare_name_paths = {
        ACES_FIT_AUTHORITY,
        SOURCE_ROOT / "app" / "config" / "config_schema.h",
        SOURCE_ROOT / "app" / "config" / "config_loader.cpp",
        SOURCE_ROOT / "app" / "config" / "session_config_adapter.cpp",
    }
    coefficient_pattern = aces_fit_sequence()
    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if LEGACY_ACES_TONEMAP_IDENTIFIER.search(code):
            errors.append(f"{relative(path)} retains the falsely named ACES tonemapper")
        if FALSE_SPECTRAL_ACES_API.search(code):
            errors.append(
                f"{relative(path)} relabels absolute spectral tristimulus values as ACES"
            )
        if DETACHED_XYZ_AP0_MATRIX.search(code):
            errors.append(
                f"{relative(path)} retains an ungoverned XYZ-to-AP0 matrix"
            )
        if path in bare_name_paths and BARE_ACES_CONFIG_LITERAL.search(document):
            errors.append(f"{relative(path)} accepts or defaults the unrepresented bare ACES name")
        if path != ACES_FIT_AUTHORITY and coefficient_pattern.search(code):
            errors.append(f"{relative(path)} reimplements the ACES-fit coefficient set")
    return errors


def verify_aces_contract_policy() -> None:
    coefficients = " ".join(value + "f" for value in ACES_FIT_COEFFICIENTS)
    valid = {
        ACES_FIT_AUTHORITY: (
            "AcesFit ParseTonemapType kDefaultTonemapperName kTonemapNames "
            "SupportedTonemapperNames "
            "std::isfinite "
            + coefficients
            + " // knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve"
        ),
        **{path: marker for path, marker in ACES_FIT_CONSUMERS.items()},
        ACES_FIT_CLI: "Tonemapper: ACESFit; not an ACES Output Transform",
        SPECTRAL_RADIANCE_FACADE: "struct SpectralRadiance {};",
    }
    if aces_contract_errors(valid):
        raise RuntimeError("ACES policy rejected the explicit represented/absent contracts")

    duplicated = dict(valid)
    consumer = next(iter(ACES_FIT_CONSUMERS))
    duplicated[consumer] += "\n" + coefficients
    if not aces_contract_errors(duplicated):
        raise RuntimeError("ACES-fit policy accepted a production coefficient copy")

    detached = dict(valid)
    detached[consumer] = "independent_tonemap();"
    if not aces_contract_errors(detached):
        raise RuntimeError("ACES-fit policy accepted a detached production consumer")

    legacy = dict(valid)
    legacy[ACES_FIT_AUTHORITY] += "\nTonemapType::Aces"
    if not aces_contract_errors(legacy):
        raise RuntimeError("ACES-fit policy accepted the falsely named internal selector")

    bare_alias = dict(valid)
    bare_alias[ACES_FIT_AUTHORITY] += '\nif (name == "ACES") return AcesFit;'
    if not aces_contract_errors(bare_alias):
        raise RuntimeError("ACES-fit policy accepted bare ACES as a represented config name")

    false_help = dict(valid)
    false_help[ACES_FIT_CLI] = "Tonemapper: ACES, Reinhard"
    if not aces_contract_errors(false_help):
        raise RuntimeError("ACES-fit policy accepted CLI advertising of an ACES transform")

    false_spectral_api = dict(valid)
    false_spectral_api[SPECTRAL_RADIANCE_FACADE] += "\nstruct Aces {}; Aces ToAces() const;"
    if not aces_contract_errors(false_spectral_api):
        raise RuntimeError("ACES policy accepted absolute spectral values labeled as ACES")

    detached_ap0 = dict(valid)
    detached_ap0[SPECTRAL_RADIANCE_FACADE] += (
        "\n1.0498110175; -0.4959030231; 1.3733130458; 0.9912520182;"
    )
    if not aces_contract_errors(detached_ap0):
        raise RuntimeError("ACES policy accepted the detached XYZ-to-AP0 matrix")


def blackbody_laws_authority_errors(documents: dict[Path, str]) -> list[str]:
    errors: list[str] = []
    authority = documents.get(BLACKBODY_LAWS_AUTHORITY)
    if authority is None:
        return ["the host blackbody-laws authority is missing"]

    authority_code = CPP_NON_CODE.sub(" ", authority)
    for marker in (
        "TryPlanckSpectralRadiancePerMetre",
        "TryWienPeakWavelength",
        "TryStefanBoltzmannExitance",
        "kPlanckC1",
        "kPlanckC2",
        "kWienB",
        "kStefanBoltzmann",
    ):
        if marker not in authority_code:
            errors.append(f"the host blackbody-laws authority omits {marker}")
    if len(PLANCK_EXPM1.findall(authority_code)) != 1:
        errors.append("the host blackbody-laws authority must own one Planck expm1")

    required_include = '#include "sirius/core/spectral/blackbody_laws.h"'
    for path in BLACKBODY_LAWS_CONSUMERS:
        document = documents.get(path)
        if document is None:
            errors.append(f"blackbody-law consumer is missing: {relative(path)}")
            continue
        if (
            required_include not in document
            or "TryPlanckSpectralRadiancePerMetre"
            not in CPP_NON_CODE.sub(" ", document)
        ):
            errors.append(
                f"{relative(path)} does not delegate to the host Planck-law authority"
            )

    constants = documents.get(PHYSICAL_CONSTANTS_AUTHORITY)
    if constants is None or len(
        WIEN_CONSTANT_LITERAL.findall(CPP_NON_CODE.sub(" ", constants))
    ) != 1:
        errors.append("the physical-constants authority must own the Wien constant once")

    for path, document in documents.items():
        code = CPP_NON_CODE.sub(" ", document)
        if path != BLACKBODY_LAWS_AUTHORITY and PLANCK_EXPM1.search(code):
            errors.append(f"{relative(path)} reimplements the host Planck law")
        if (
            path not in {BLACKBODY_LAWS_AUTHORITY, PHYSICAL_CONSTANTS_AUTHORITY}
            and PLANCK_CONSTANT_IDENTIFIER.search(code)
        ):
            errors.append(f"{relative(path)} consumes Planck constants outside their authority")
        if path != PHYSICAL_CONSTANTS_AUTHORITY and WIEN_CONSTANT_LITERAL.search(code):
            errors.append(f"{relative(path)} copies the Wien displacement constant")
        if "StefanBoltzmannRadiance" in code:
            errors.append(f"{relative(path)} retains the mislabeled Stefan-Boltzmann quantity")
    return errors


def verify_blackbody_laws_authority_policy() -> None:
    valid = {
        BLACKBODY_LAWS_AUTHORITY: (
            "TryPlanckSpectralRadiancePerMetre std::expm1(exponent) kPlanckC1 kPlanckC2 "
            "TryWienPeakWavelength kWienB TryStefanBoltzmannExitance kStefanBoltzmann"
        ),
        **{
            path: (
                '#include "sirius/core/spectral/blackbody_laws.h"\n'
                "TryPlanckSpectralRadiancePerMetre"
            )
            for path in BLACKBODY_LAWS_CONSUMERS
        },
        PHYSICAL_CONSTANTS_AUTHORITY: "kWienB = 2.897771955e-3",
    }
    if blackbody_laws_authority_errors(valid):
        raise RuntimeError("blackbody-laws policy rejected the single host authority")

    consumer = BLACKBODY_LAWS_CONSUMERS[0]
    duplicated = dict(valid)
    duplicated[consumer] += "\nstd::expm1(exponent)"
    if not blackbody_laws_authority_errors(duplicated):
        raise RuntimeError("blackbody-laws policy accepted a production Planck duplicate")

    reconstructed = dict(valid)
    reconstructed[consumer] += "\nkPlanck * kSpeedOfLight / kBoltzmann"
    if not blackbody_laws_authority_errors(reconstructed):
        raise RuntimeError("blackbody-laws policy accepted reconstructed Planck constants")

    detached = dict(valid)
    detached[consumer] = "independent_planck_law();"
    if not blackbody_laws_authority_errors(detached):
        raise RuntimeError("blackbody-laws policy accepted a detached production consumer")

    copied_wien = dict(valid)
    copied_wien[consumer] += "\n2.897771955e-3"
    if not blackbody_laws_authority_errors(copied_wien):
        raise RuntimeError("blackbody-laws policy accepted a copied Wien constant")

    mislabeled = dict(valid)
    mislabeled[consumer] += "\nStefanBoltzmannRadiance(temperature)"
    if not blackbody_laws_authority_errors(mislabeled):
        raise RuntimeError("blackbody-laws policy accepted a mislabeled physical quantity")
