#!/usr/bin/env python3
"""Freeze original model-selected CPU camera controls without execution."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PRODUCT = "4bcd9875cd9308deabed06a5fa2f17d2321b0cbc"
RUNNER = "a1d6e75299a2bf250e5ce7d44ad5965f8acb4299"
OWNER_ARCHIVE = ROOT / "attestations/scientific-core/4bcd987/p6-cpu-writer-separated-controls"
OWNER_MANIFEST = "0124c6fc757670470ccfa1909843d8de96beebd25a910fdf4a48b0b981b3add3"
OWNER_RAW = OWNER_ARCHIVE / "files/out/p6-cpu-writer-separated-controls"
CORE = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
CORE_SHA = "c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58"
RENDER = ROOT / "bin/linux-gcc/src/sirius/app/sirius_render_tests"
RENDER_SHA = "2489bc56bd0c1c53ab01da021119d6292d8fdb6f51037b82819d3dcdf6e87ab4"
CASE = "original-24-cpu-camera-controls"
EXTRA = "CameraContinuousProjection.AngularProjectionPreservesGeometricZerosForDetectorOffsets"
PREFIXES = ("CameraContinuousProjection.", "CameraFilmDifferential.", "CameraSampling.",
            "CameraWorldlineTest.", "ThinLensCameraTest.", "ObserverFrameTests.")
SOURCES = [
    "tests/core/camera_continuous_projection_test.cpp",
    "tests/core/camera_film_differential_test.cpp",
    "tests/core/camera_sampling_test.cpp",
    "tests/core/camera_test.cpp",
    "tests/core/geodesic_live_path_test.cpp",
    "tests/render/camera_aberration_test.cpp",
]
HEADERS = ["src/sirius/core/" + name for name in (
    "camera.h", "camera_launch.h", "camera_sampling.h", "observer_frame.h",
    "celestial_tangent_basis.h", "metrics/kerr_schild_family.h")]


def identity(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def key(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def write(name, value):
    (WORK / name).write_text(json.dumps(value, indent=2) + "\n")


assert not (WORK / "preparation.json").exists() and not (WORK / CASE).exists()
assert git("rev-parse", "HEAD").decode().strip() == RUNNER and not git("status", "--porcelain").strip()
changed = git("diff", "--name-only", PRODUCT, RUNNER).decode().splitlines()
assert set(changed) == {".github/workflows/ci.yml", "scripts/diagnose-macos-retained.py",
                        "scripts/governance/build_policy.py", "scripts/macos_metal_probe.py"}
assert identity(OWNER_ARCHIVE / "manifest.json")["sha256"] == OWNER_MANIFEST
manifest = json.loads((OWNER_ARCHIVE / "manifest.json").read_text())
records = {row["path"]: {field: row[field] for field in ("bytes", "sha256")} for row in manifest["files"]}
base = json.loads((OWNER_RAW / "preparation.json").read_text())
assert identity(OWNER_RAW / "preparation.json") == records[key(OWNER_RAW / "preparation.json").split(OWNER_ARCHIVE.name + "/", 1)[1]]
model = json.loads((ROOT / "tests/operating_model.json").read_text())
criteria = {row["id"]: row for row in model["acceptance_criteria"]}
p5 = [row["name"] for row in criteria["P5"]["evidence"]
      if row["kind"] == "gtest" and row["name"].startswith(PREFIXES)]
assert len(p5) == len(set(p5)) == 23 and EXTRA not in p5
for criterion in ("P2", "P3"):
    assert {"kind": "gtest", "name": EXTRA} in criteria[criterion]["evidence"]
names = sorted([*p5, EXTRA])
wrapper_ast = ast.parse((WORK / "run_camera.py").read_text())
wrapper_names = next(ast.literal_eval(node.value) for node in wrapper_ast.body
                     if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "NAMES" for t in node.targets))
assert set(names) == wrapper_names and len(names) == 24
pins = {name: pin for name, pin in base["pins"].items()
        if name.startswith(("bin/", "/usr/", "attestations/alignment/", "attestations/fresh-session-review/"))}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
core_recovery = ROOT / ("attestations/alignment/4bcd987/gcc-revision-synchronisation/prior-versions/" + CORE_SHA)
render_recovery = ROOT / ("attestations/alignment/4bcd987/gcc-revision-synchronisation/current-versions/" + RENDER_SHA)
assert identity(CORE) == identity(core_recovery) == {"bytes": 5199120, "sha256": CORE_SHA}
assert identity(RENDER) == identity(render_recovery) == {"bytes": 12208240, "sha256": RENDER_SHA}
accepted_path = ROOT / "attestations/alignment/4bcd987/gcc-revision-synchronisation/files/out/gcc-revision-synchronisation/accepted-result.json"
accepted = json.loads(accepted_path.read_text())
assert accepted["passed"] and accepted["source_revision"] == PRODUCT
for name in ("run_owned.py", "owned_guard.py"):
    original = OWNER_RAW / name
    assert identity(original) == records[str(original.relative_to(OWNER_ARCHIVE))]
    shutil.copy2(original, WORK / name)
source_records = {}
for name in [*SOURCES, *HEADERS, "tests/operating_model.json", "tests/core/CMakeLists.txt", "tests/render/CMakeLists.txt", "tests/labels/CTestLabels.cmake"]:
    path = ROOT / name
    original = git("show", PRODUCT + ":" + name)
    assert original == path.read_bytes(), name
    source_records[name] = {**identity(path), "git_blob": git("rev-parse", PRODUCT + ":" + name).decode().strip()}
    target = WORK / "reference" / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[key(path)] = identity(path)
    pins[key(target)] = identity(target)
compile_commands = json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
compiles = []
source_for_name = {}
for source in SOURCES:
    rows = [row for row in compile_commands if Path(row["file"]).resolve() == ROOT / source]
    assert len(rows) == 1, source
    row = rows[0]
    assert row["directory"] == str(ROOT / "bin/linux-gcc")
    for flag in ("/usr/bin/g++-14", "-DSIRIUS_CONTRACT_MODE=2", "-O3", "-DNDEBUG", "-std=c++2c", "-Werror", "tests/support/fail_on_skip.h"):
        assert flag in row["command"], (source, flag)
    target = "sirius_render_tests" if source.startswith("tests/render/") else "sirius_core_tests"
    assert target + ".dir/" in row["output"] and row["output"].endswith(Path(source).name + ".o")
    compiles.append(row)
    for suite, case in re.findall(r"TEST\(\s*(\w+)\s*,\s*(\w+)\s*\)", (ROOT / source).read_text()):
        full = suite + "." + case
        if full in names:
            assert full not in source_for_name
            source_for_name[full] = source
assert set(source_for_name) == set(names)
registrations = []
registered = {}
for relative, producer, expected_count in (
    ("bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake", CORE, 23),
    ("bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake", RENDER, 1),
):
    lines = (ROOT / relative).read_text().splitlines()
    chosen = []
    for name in names:
        selected = [line for line in lines if line.startswith("add_test([=[" + name + "]=]")]
        if not selected:
            continue
        assert len(selected) == 1 and str(producer) + " " in selected[0]
        assert "--gtest_filter=" + name + "]==]" in selected[0]
        props = [line for line in lines if line.startswith("set_tests_properties([=[" + name + "]=]")]
        assert len(props) == 1 and "LABELS Correctness" in props[0] and "TIMEOUT" not in props[0]
        assert name not in registered
        record = {"name": name, "source": source_for_name[name], "registration_file": relative,
                  "add_test": selected[0], "properties": props[0], "producer": key(producer)}
        chosen.append(record)
        registered[name] = record
    assert len(chosen) == expected_count
    registrations.extend(chosen)
assert set(registered) == set(names)
ninja = (ROOT / "bin/linux-gcc/build.ninja").read_text().splitlines()
links = [line for line in ninja if line.startswith(("build tests/core/sirius_core_tests ", "build src/sirius/app/sirius_render_tests "))
         and ": CXX_EXECUTABLE_LINKER__" in line]
assert len(links) == 2
for compile in compiles:
    assert sum(compile["output"] in line for line in links) == 1

# These descriptions retain the original finite domains and assertions; source
# snapshots supply every unchanged detail rather than a new executable oracle.
scopes = {
    "CameraContinuousProjection.GeometryAndDifferentialAcrossOutputCrop": "1280x640; yaw .13/pitch -.21/roll .09; all three lenses; four film points (-12.25,-8.75),(1291.3,649.1),(437.375,255.125),(640,320); pupil(.2,.7), thin focus30. Independent normalized plane-rotation/fixed-focus/spherical geometry 5e-15; central FD h=.002 at 3e-12; cross-product density 2e-20.",
    "CameraContinuousProjection.FixedPupilAndObserverArePreserved": "ThinLens1280x640 focus30 beta(.2,-.1,.3), four original samples, film(-3.25,650.75). Exact pupil/origin/beta preservation; fixed-focus plane intersections 2e-13.",
    "CameraContinuousProjection.FisheyeCropMaskCentreAndAntipodalRim": "100x100/FOV60 at (150,50),(350,50),(350.1,50), then FOV360 rim(100,50)/centre(50,50). Original active/represented/differential-mask checks; centre direction exact -1 and density pi_float^2*.02^2 at 2e-18.",
    "CameraContinuousProjection.RefinementDeclinesLostCoordinatesAndDirections": "Default Pinhole offsets(1e-8,-1e-8) at(960,540) remain distinct; (1e-20,0),(1e-20,.25) and denorm_min at(0,540) decline as Unrepresentable.",
    "CameraContinuousProjection.CropTranslationDoesNotRenormalizeOrMaskProjection": "Three lenses,1280x640 versus1600x640; x=-16,.25,1279.75,1296 translated+160 at y170.25,pupil(.2,.7). Direction5e-16 and density2e-20.",
    "CameraContinuousProjection.OriginalPacketsKeepLegacyValuesWithExplicitSmoothRoundoff": "1920x1080; yaw .17/pitch-.11/roll.07,beta_x.2; all three lenses/four original samples/three pixels(0,0),(1919,1079),(713,431), thin focus30. Legacy before/after ray exact, pupil exact; maximum legacy/smooth direction difference >1e-10 and <1e-6 (compatibility budget only).",
    "CameraContinuousProjection.InvalidAndUnsupportedProjectionPublishesNoRay": "Default Pinhole: NaN/+Inf film, pupil u=-1/1/NaN ->InvalidInput; max_double film->Arithmetic; original custom camera->Unsupported, no ray.",
    EXTRA: "ThinLens4x2/FOV2/focus50, film(.5,.5),pupil(.2,.7). Exact angular Jacobian p[1][1]=0 and nonzero determinant; invert original sigma=(2*pi/180)/2*.3 offsets for z=(4,0),(0,4), both projections represented.",
    "CameraFilmDifferential.IndependentGeometryAndActualRayStencils": "Three lenses1200x800 yaw.31/pitch-.27/roll.19; pixels(121,231),(731,437),(990,681),sample(.37,.63,.21,.82). Independent double normalized geometry FD h=.01:2e-9+2e-6*abs(derivative); actual public float-ray stencil offsets[-4,-2,2,4]/24:8e-8. Celestial reconstruction/longitudinal orthogonality2e-17; pupil exact; positive absolute/signed density.",
    "CameraFilmDifferential.FullLaunchDifferentiatesDisplacedMovingKerrObserver": "Kerr mass1/spins[-.9,0,.9]; observer r12/theta1.1/phi.7,beta(.31,-.17,.21),ThinLens1920x1080/focus20,film(704.25,367.75),pupil(.73,.29). Four independent columns film x/y and Cartesian pupil right/up; complete launches at offsets[-2h,-h,h,2h],h=.002, with analytic phase_space reset and fixed focus recovered for pupil offsets. Five-point displacement/coordinate tangent derivatives2e-10; central null/frequency2e-14; differentiated null2e-13; non-frozen frame contribution>1e-4. Reference still uses original LaunchCameraRay for varied rays, not an independent launch implementation.",
    "CameraFilmDifferential.FilmColumnsRecoverAngularScreenAndPinholeHasNoPupilColumns": "Kerr mass1/spin.6; default Pinhole with beta(.2,-.1,.3),yaw.3/pitch.2/roll-.4,film(437.25,720.125). Angular columns match observer screen2e-14; angular displacements and both pupil columns exact zero. NaN pupil derivative declines launch; singular film_to_angle declines angular inversion.",
    "CameraFilmDifferential.PerspectiveAndFixedPupilSolidAngleLaws": "Default Pinhole/ThinLens; pixel(283,729),film sample(.29,.71),pupils(.5,0),(.19,.93). Independent plane-area law dOmega=normal_distance*dA/distance^3 at relative3e-6, positive sign; centre-pupil equivalence at(127,349) relative3e-6.",
    "CameraFilmDifferential.EquidistantCentreLimitAndAntipodalMask": "Fisheye1000x1000/FOV360,pixels(500,500),(501,500),(672,793),zero sample. Independent equidistant a*sin(a*r)/r density (centre a^2) relative3e-6; centre derivatives1e-15 and radial exact zeros. Original antipodal active/no differential and outside-optics inactive/no differential checks.",
    "CameraFilmDifferential.PixelUnitsAndIntegratedRectangleSolidAngle": "Pinhole320x240 composite midpoint sum over every original pixel versus rectangle 4*atan(a*b/sqrt(1+a^2+b^2)) at2e-5. Separate binary-exact512x256->1024x512 samples(123,77,.25,.25)->(246,154,.5,.5): EXPECT_DOUBLE_EQ density*4 and derivatives*2.",
    "CameraFilmDifferential.PreservesNominalRaysAndExistingOrientationSemantics": "Default three lenses at(297,381),four original samples; rotated yaw.31/pitch-.27/roll.19. Before/after rays exact; Pinhole direction changes/density relative3e-6; original ThinLens/Fisheye orientation semantics preserve rays/Jacobians exactly.",
    "CameraFilmDifferential.GovernedFourSampleLaunchAndSingleObserverBoost": "ThinLens r30/theta80deg/focus30,beta(.1,.02,-.01),three pixels(169,229),(1167,586),(874,555),four original samples =>12 launches. Same CheckDifferential bounds as geometry/stencil case; stationary and moving rest-camera maps exact. Observer-screen/rest derivative reconstruction2e-17 and independently differenced boosted geometric null launches h=.01 at2e-9+2e-6*abs(expected).",
    "CameraFilmDifferential.CustomCameraDeclinesUnavailableDifferential": "Original custom Pinhole-delegating camera at(50,60): nominal observer ray active, differential unavailable; no fabricated Jacobian required.",
    "ObserverFrameTests.MovingCameraLaunchMatchesIndependentLorentzTransform": "Minkowski Cartesian unit seeds, beta(.25,-.4,.1),screen n=(-.8,.3,sqrt(.27)). Independent Lorentz q'=(1,-n), gamma and beta-dot-n equations compare all four photon components2e-14.",
    "ObserverFrameTests.KerrCameraRayIsPastNullAndLaunchFrequencyIsOne": "Kerr mass1/spin.9,event Cartesian(30,2,4),beta(.7,-.2,.1),screen(-.9,.3,sqrt(.1)). Past-directed time<0; null,unit launch frequency,tetrad signature,screen orthonormality/observer/null orthogonality all3e-13.",
    "CameraWorldlineTest.RestScreenRayAndWorldlineComposeOverLensModels": "64x64,beta(.3,-.2,.1),all three lenses,pixel(20,24),sample(.5,.5). Bound/rest direction components exact and beta fields EXPECT_DOUBLE_EQ. Camera construction/projection only; no RenderSession or frame rendering.",
    "ThinLensCameraTest.CentreApertureSharesPinholeProjectionAndFieldOfView": "200x100/FOV60/focal50/aperture2.8/focus30 versus Pinhole default focus,GenerateRay(150,50,.25,0). Centred-pupil directions1e-6; FOV30 must narrow azimuth magnitude.",
    "ThinLensCameraTest.PupilSamplesSpanFiniteLaunchAndConvergeAtFocus": "320x180/FOV52/focal50/aperture2/focus30,pixel(237,71),film(.35,.65),pupils(.25,.64),(.70,.36). Pupil radius>0 and<=original aperture radius; forward direction>0; independent focus-plane intersections5e-6.",
    "ThinLensCameraTest.OversizedPupilDeclinesTheLocalTangentPlane": "r1/focus1/focal100/aperture1 violates original r_pupil<=.1*min(r_observer,d_focus). Configuration issue and original precondition death assertion; no new boundary sample.",
    "CameraSampling.SequenceDoesNotCollapseFilmAndPupilDimensions": "Exactly4096 emitted/observed original deterministic samples; every film/pupil coordinate strictly(0,1). Mean[(film_u-pupil_u)^2+(film_v-pupil_v)^2] within.01 of1/3; no random-uniform/all-sequence claim.",
}
assert set(scopes) == set(names)
local = [WORK / name for name in ("prepare.py", "run_camera.py", "run_owned.py", "owned_guard.py")]
for path in local:
    ast.parse(path.read_text(), filename=str(path))
extra_paths = [CORE, core_recovery, render_recovery, OWNER_ARCHIVE / "manifest.json", OWNER_ARCHIVE / "disposition.json",
               OWNER_RAW / "run_owned.py", OWNER_RAW / "owned_guard.py", *local]
for path in extra_paths:
    pins[key(path)] = identity(path)
ctest = re.search(r"^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$", (ROOT / "bin/linux-gcc/CMakeCache.txt").read_text(), re.M)
assert ctest
actions = {CASE: {"revision": RUNNER, "kind": "check", "expected_exit": 0,
    "limits": {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152},
    "scope": "Twenty-four unchanged original finite CPU camera controls only; no device dispatch, frame, full P5/P2/P3, Mandatory or release qualification.",
    "argv": [ctest.group(1), "--test-dir", "bin/linux-gcc", "-V", "--output-on-failure", "--no-tests=error",
             "--output-junit", "{folder}/ctest.xml", "--output-log", "{folder}/ctest.log", "-R",
             "^(" + "|".join(re.escape(name) for name in names) + ")$"],
    "frozen_inputs": list(pins) + [key(WORK / "actions.json"), key(WORK / "preparation.json")]}}
write("actions.json", actions)
pins[key(WORK / "actions.json")] = identity(WORK / "actions.json")
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": PRODUCT, "runner_source_revision": RUNNER, "diagnostic_only_difference": changed,
    "executed": False, "expected_names": names, "scope": actions[CASE]["scope"],
    "ledger_scope": {"P5_original_selected_names": p5, "P2_P3_extra": EXTRA, "unique_count": 24,
                     "not_claimed": ["all finite timelike velocities", "all lens/film/pupil domains", "independent Slang projection", "retained GPU parity", "5616x4096 frame", "full P5/P2/P3", "Mandatory", "release"]},
    "finite_acceptance_scope": scopes,
    "reference": {"kind": "Existing original geometric, finite-difference, plane-solid-angle, rectangle, Lorentz and sequence assertions; no new oracle or stronger threshold.",
        "equations": {"normalized_direction": "n=q/|q|; its tangent is (I-n*n^T)*dq/|q|; existing independent test geometry is normalized in double then finite-differenced, not a copied production derivative.",
                      "fixed_pupil": "Intersect the original fixed pupil ray with its requested focus plane.",
                      "solid_angle": "dOmega=normal_distance*dA/distance^3; rectangle=4*atan(a*b/sqrt(1+a^2+b^2)).",
                      "moving_frame": "Original independent Lorentz transform of physical photon q'=(1,-n).",
                      "four_columns": "Original five-point complete-launch variations in film x/y pixels and physical Cartesian pupil right/up; stochastic disk-map differentiation is outside this test."},
        "original_source_snapshots": source_records,
        "default_camera_config": {"t": 0, "r": 50, "theta": "pi/2", "phi": 0, "yaw_pitch_roll": [0,0,0],
            "beta_xyz": [0,0,0], "fov": 60, "focal_length": 50, "aperture": 2.8, "focus_distance": 50, "width": 1920, "height": 1080},
        "representation_boundary": "Original float lens coefficients and legacy-ray roundoff budgets are retained; these finite checks do not establish all-configuration bounds."},
    "source_and_registration_join": {"test_sources": SOURCES, "production_headers": HEADERS,
        "core_producer": identity(CORE), "core_recovery_archive": key(core_recovery),
        "render_producer": identity(RENDER), "render_recovery_archive": key(render_recovery),
        "accepted_source_build_evidence": key(accepted_path), "runtime_resources_bound": base["source_and_registration_join"]["runtime_resource_files"],
        "cache_required": base["source_and_registration_join"]["cache_required"], "original_compile_commands": compiles,
        "original_link_ninja_lines": links, "registrations": registrations,
        "selected_core_cases": 23, "selected_render_cases": 1, "original_test_TIMEOUT_added": False},
    "owner_reuse": {"archive": key(OWNER_ARCHIVE), "manifest_sha256": OWNER_MANIFEST,
                    **{name: {"source": key(OWNER_RAW / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")}},
    "private_tmpdir": {"path": str(WORK / "temporary"), "created_only_at_execution": True, "mode": "0700", "owner_uid": os.getuid(), "preserve_natural_residuals": True},
    "pins": pins,
})
print(json.dumps({"prepared_without_execution": True, "names": len(names), "core": 23, "render": 1, "pins": len(pins),
                  "preparation": identity(WORK / "preparation.json"), "command": ["python3", "-B", key(WORK / "run_camera.py")]}, indent=2))
