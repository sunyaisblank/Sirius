#!/usr/bin/env python3
"""Freeze nine approved original CPU radiation controls without execution."""
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
CASE = "original-nine-cpu-radiation-controls"
OWNER = ROOT / "attestations/scientific-core/4bcd987/p5-cpu-camera-current-source"
OWNER_MANIFEST = "3dce840b04e1909f9d83d7d7e8bf360ac61fc7bcf3d6ca1dbf8d898fee591b08"
OWNER_RAW = OWNER / "files/out/p5-cpu-camera-current-source"
SOURCES = ["tests/core/blackbody_planck_test.cpp", "tests/core/spectral_radiance_test.cpp",
           "tests/render/spectral_emission_test.cpp"]
HEADERS = ["src/sirius/core/" + name for name in (
    "constants.h", "spectral/blackbody_laws.h", "spectral/blackbody.h",
    "spectral/spectral_radiance.h", "spectral/colour_modes.h",
    "relativistic_transfer.h", "kerr_orbits.h")]


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
assert identity(OWNER / "manifest.json")["sha256"] == OWNER_MANIFEST
manifest = json.loads((OWNER / "manifest.json").read_text())
records = {row["path"]: {field: row[field] for field in ("bytes", "sha256")} for row in manifest["files"]}
base = json.loads((OWNER_RAW / "preparation.json").read_text())
assert identity(OWNER_RAW / "preparation.json") == records[str((OWNER_RAW / "preparation.json").relative_to(OWNER))]
wrapper = ast.parse((WORK / "run_radiation.py").read_text())
names = sorted(next(ast.literal_eval(node.value) for node in wrapper.body
                    if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "NAMES" for t in node.targets)))
assert len(names) == len(set(names)) == 9
model = json.loads((ROOT / "tests/operating_model.json").read_text())
p4 = next(row for row in model["acceptance_criteria"] if row["id"] == "P4")
assert set(names) <= {row["name"] for row in p4["evidence"] if row["kind"] == "gtest"}
pins = {name: pin for name, pin in base["pins"].items()
        if name.startswith(("bin/", "/usr/", "attestations/alignment/", "attestations/fresh-session-review/"))}
for name, pin in pins.items():
    assert identity(ROOT / name) == pin, name
joins = base["source_and_registration_join"]
core = ROOT / "bin/linux-gcc/tests/core/sirius_core_tests"
render = ROOT / "bin/linux-gcc/src/sirius/app/sirius_render_tests"
assert identity(core) == joins["core_producer"] == identity(ROOT / joins["core_recovery_archive"])
assert identity(render) == joins["render_producer"] == identity(ROOT / joins["render_recovery_archive"])
assert joins["core_producer"] == {"bytes": 5199120, "sha256": "c0c076aac164a3b8c917bd4c66e7d388a4a9f70a2f8a56b40e2a088429886f58"}
assert joins["render_producer"] == {"bytes": 12208240, "sha256": "2489bc56bd0c1c53ab01da021119d6292d8fdb6f51037b82819d3dcdf6e87ab4"}
for name in ("run_owned.py", "owned_guard.py"):
    original = OWNER_RAW / name
    assert identity(original) == records[str(original.relative_to(OWNER))]
    shutil.copy2(original, WORK / name)
snapshots = {}
for name in [*SOURCES, *HEADERS, "tests/operating_model.json", "tests/core/CMakeLists.txt", "tests/render/CMakeLists.txt", "tests/labels/CTestLabels.cmake"]:
    path = ROOT / name
    assert git("show", PRODUCT + ":" + name) == path.read_bytes(), name
    snapshots[name] = {**identity(path), "git_blob": git("rev-parse", PRODUCT + ":" + name).decode().strip()}
    target = WORK / "reference" / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    pins[name] = identity(path)
    pins[key(target)] = identity(target)
compiles = []
source_for_name = {}
commands = json.loads((ROOT / "bin/linux-gcc/compile_commands.json").read_text())
for source in SOURCES:
    rows = [row for row in commands if Path(row["file"]).resolve() == ROOT / source]
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
for relative, producer, count in (
    ("bin/linux-gcc/tests/core/sirius_core_tests[1]_tests.cmake", core, 8),
    ("bin/linux-gcc/tests/render/sirius_render_tests[1]_tests.cmake", render, 1),
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
        chosen.append({"name": name, "source": source_for_name[name], "registration_file": relative,
                       "add_test": selected[0], "properties": props[0], "producer": key(producer)})
    assert len(chosen) == count
    registrations.extend(chosen)
assert {row["name"] for row in registrations} == set(names)
links = [line for line in (ROOT / "bin/linux-gcc/build.ninja").read_text().splitlines()
         if line.startswith(("build tests/core/sirius_core_tests ", "build src/sirius/app/sirius_render_tests "))
         and ": CXX_EXECUTABLE_LINKER__" in line]
assert len(links) == 2 and all(sum(row["output"] in line for line in links) == 1 for row in compiles)
nist_path = WORK / "reference/nist-2022-allascii.txt"
nist_provenance_path = WORK / "reference/nist-2022-provenance.json"
nist = json.loads(nist_provenance_path.read_text())
assert identity(nist_path) == {field: nist[field] for field in ("bytes", "sha256")}
assert nist["url"] == "https://physics.nist.gov/cuu/Constants/Table/allascii.txt" and nist["edition"] == "2022 CODATA adjustment"
assert "2022 CODATA adjustment" in nist_path.read_text()
for row in nist["primary_rows"].values():
    assert row in nist_path.read_text().splitlines()
values = {name: row[60:85].replace(" ", "") for name, row in nist["primary_rows"].items()}
assert values["Planck constant"] == "6.62607015e-34"
assert values["speed of light in vacuum"] == "299792458"
assert values["Boltzmann constant"] == "1.380649e-23"
assert "5.670374419...e-8" == values["Stefan-Boltzmann constant"]
assert "2.897771955...e-3" == values["Wien wavelength displacement law constant"]

scopes = {
    "SpectralValidationTests.PlanckSpectralRadianceMatchesIndependentCodataEquation": {
        "reference": "Independent 2*h*c^2/[lambda^5*expm1(h*c/(lambda*k*T))] with original literal SI h,c,k, not production Planck coefficients.",
        "temperatures_K": [1200, 3000, 5778, 10000, 1e7], "wavelengths_m": [100e-9, 380e-9, 550e-9, 780e-9, 10e-6],
        "sample_pairs": 25, "assertion": "Original relative error <2e-14 (denominator max(abs(expected),1e-300)); no all-wavelength/all-temperature bound."},
    "SpectralValidationTests.WienPeakMatchesIndependentNumericalPlanckMaximum": {
        "temperatures_K": [1200, 3000, 5778, 10000], "reference": "Original 120-iteration golden-section maximization of independent Planck on [50e-9,20e-6] m.",
        "assertion": "Original relative error <2e-8; no new root solver or executable reference."},
    "SpectralValidationTests.StefanBoltzmannExitanceMatchesIndependentHemisphericPlanckIntegral": {
        "temperatures_K": [1200, 3000, 5778, 10000], "reference": "Original 32768-panel midpoint integral in log(lambda), log limits -20.72326583694641 to -6.907755278982137 (1 nm to 1 mm), independent Planck times lambda and hemisphere factor pi.",
        "assertion": "Original relative error <2e-7 against production sigma*T^4; finite integration bounds retained."},
    "SpectralValidationTests.BlackbodyLawsRejectUnrepresentedDomains": {
        "invalid_values": [0, -1, "+infinity", "NaN"], "reference": "Original invalid wavelength/temperature refusals; deep Wien lambda1e-20 m/T1 K exact zero; finite extreme Rayleigh-Jeans lambda1e100 m/T1e300 K compares original c1/c2*1e-100 at relative<2e-13. This extreme-limit expected coefficient ratio is production-derived.",
        "assertion": "Original overflow refusal at lambda1e-300 m,T=c2/1e-297 and Stefan(max_double); no new domain points."},
    "SpectralValidationTests.KerrDiskTransferMatchesKillingFieldContraction": {
        "domain": {"mass": 1, "signed_spin": 0.7, "equatorial_radius": 6, "E": 1, "Lz": 2, "observer_frequency": 0.93},
        "reference": "Original separately written equatorial BL g_tt/g_tphi/g_phiphi, circular Omega, normalized u^t and nu=u^t*(E-Omega*Lz); same independent contraction for ZAMO Omega=-g_tphi/g_phiphi.",
        "assertion": "Original absolute 2e-15 emitter/zamo frequencies and full_g/zamo_g; one fixed tuple, no arbitrary-latitude claim."},
    "SpectralValidationTests.ZamoBranchRemainsTimelikeInsideTheErgosphere": {
        "domain": {"mass": 1, "signed_spin": 0.998, "equatorial_radius": 1.5, "E": 1, "Lz": 1.5, "observer_frequency": 1},
        "assertion": "Original static ratio(gtt_emit=1/3,gtt_observer=-1) declines; Kerr transfer exists, zamo_frequency>0 and zamo_g finite. No additional numeric oracle inferred."},
    "SpectralValidationTests.OpticallyThinGreyLayerPreservesFirstOrderEmission": {
        "domain": {"delta_tau": 1e-14, "observed_source": [1, 0.5, 0.25], "max_tau": 10},
        "reference": "Original first-order 1-exp(-tau)=tau+O(tau^2) ratios at the one thin layer; no E2/polarised-transport rerun.",
        "assertion": "Original positive absorbed fraction/emission, absorbed/tau near1 at1e-14, negative tau declines; accumulated RGB/tau near[1,.5,.25] at1e-14."},
    "SpectralRadianceTest.BlackbodyBinsDelegateToPlanckAuthorityAndRejectInvalidTemperature": {
        "domain": {"temperature_K": 5778, "bins": 32, "bin_centres_nm": "380+(bin+.5)*12.5 for bin0..31"},
        "reference": "Facade/delegation and unit conversion against production TryPlanckSpectralRadiancePerMetre, not an independent 32-bin physical oracle.",
        "assertion": "Original EXPECT_DOUBLE_EQ L[bin]=production_per_metre*1e-9 (per nm); temperatures0,-1,+Inf,NaN yield exact total energy0."},
    "SpectralEmissionTest.BolometricDiskAuthorityAppliesExactlyOneGFourthFactor": {
        "domain": {"emitted_float": 0.37, "g_float": 1.25}, "reference": "Original algebraic bolometric authority check: invariant I_nu/nu^3 contributes g^3 and dnu_observer=g*dnu_emitter contributes one g, hence g^4.",
        "assertion": "Original EXPECT_FLOAT_EQ ObservedBolometricIntensity(emitted,g)=emitted*g*g*g*g at one tuple; no integrated disk/frame or pixel accuracy claim."},
}
assert set(scopes) == set(names)
historical = ROOT / "attestations/scientific-core/5e70607/final-build-identity-5e70607.json"
prior = json.loads(historical.read_text())
assert prior["artifacts"][key(core)] != identity(core)["sha256"]
assert prior["artifacts"][key(render)] != identity(render)["sha256"]
local = [WORK / name for name in ("prepare.py", "run_radiation.py", "run_owned.py", "owned_guard.py")]
for path in local:
    ast.parse(path.read_text(), filename=str(path))
extra = [*local, nist_path, nist_provenance_path, OWNER / "manifest.json", OWNER / "disposition.json",
         OWNER_RAW / "preparation.json", OWNER_RAW / "run_owned.py", OWNER_RAW / "owned_guard.py", historical,
         historical.parent / "final-integrated-cpu.xml"]
for path in extra:
    pins[key(path)] = identity(path)
ctest = re.search(r"^CMAKE_CTEST_COMMAND:INTERNAL=(.+)$", (ROOT / "bin/linux-gcc/CMakeCache.txt").read_text(), re.M)
assert ctest
action = {"revision": RUNNER, "kind": "check", "expected_exit": 0,
    "limits": {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152},
    "scope": "Nine unchanged original finite CPU radiation/transfer controls only; no device/frame, full P4/science, Mandatory or release qualification.",
    "argv": [ctest.group(1), "--test-dir", "bin/linux-gcc", "-V", "--output-on-failure", "--no-tests=error",
             "--output-junit", "{folder}/ctest.xml", "--output-log", "{folder}/ctest.log", "-R",
             "^(" + "|".join(re.escape(name) for name in names) + ")$"],
    "frozen_inputs": list(pins) + [key(WORK / "actions.json"), key(WORK / "preparation.json")]}
write("actions.json", {CASE: action})
pins[key(WORK / "actions.json")] = identity(WORK / "actions.json")
write("preparation.json", {
    "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "product_source_revision": PRODUCT, "runner_source_revision": RUNNER, "diagnostic_only_difference": changed,
    "executed": False, "expected_names": names, "scope": action["scope"], "finite_acceptance_scope": scopes,
    "ledger_scope": {"criterion": "P4", "selected_original_names": names, "selected_count": 9,
                     "excluded": "Previously completed E2/polarised witnesses, unselected P4 controls, all device and frame work."},
    "reference": {"original_source_snapshots": snapshots, "nist_2022": nist,
        "nist_scope": "Official table binds coefficient provenance only; no pixel, all-radiation, or whole-domain accuracy claim.",
        "derivations_for_review_only": {"Planck": "B_lambda=2*h*c^2/[lambda^5*(exp(h*c/(lambda*k*T))-1)] per metre and per sr.",
            "Wien": "Positive root x of5*(1-exp(-x))=x; b=(h*c/k)/x. Original executable reference remains golden-section Planck maximum.",
            "hemisphere_exitance": "pi*integral_0^infinity B_lambda dlambda=[2*pi^5*k^4/(15*h^3*c^2)]*T^4. Original executable reference remains bounded log-lambda quadrature.",
            "Kerr_frequency": "nu=-u dot p=u^t*(E-Omega*Lz), g=nu_observer/nu_emitter.",
            "bolometric": "Invariant I_nu/nu^3 plus frequency measure dnu gives g^4; one algebraic float test only."},
        "source_composition_boundary": "Root reviewed blackbody law, constants, transfer composition and facade units; finite source algebra review is not exhaustive or a machine proof. No new numerical reference implementation or threshold.",
        "historical_execution": {"xml": key(historical.parent / "final-integrated-cpu.xml"), "identity": key(historical),
            "revision": prior["revision"], "producer_hashes": prior["artifacts"],
            "reuse_boundary": "Seven selected SpectralValidation cases have historical passes at5e70607 with different ELF hashes; no matching current-ELF selected execution found in preserved XML. No relabelling as current."}},
    "source_and_registration_join": {"test_sources": SOURCES, "production_headers": HEADERS,
        "core_producer": identity(core), "core_recovery_archive": joins["core_recovery_archive"],
        "render_producer": identity(render), "render_recovery_archive": joins["render_recovery_archive"],
        "accepted_source_build_evidence": joins["accepted_source_build_evidence"], "cache_required": joins["cache_required"],
        "runtime_resources_bound": joins["runtime_resources_bound"], "original_compile_commands": compiles,
        "original_link_ninja_lines": links, "registrations": registrations,
        "selected_core_cases": 8, "selected_render_cases": 1, "original_test_TIMEOUT_added": False},
    "owner_reuse": {"archive": key(OWNER), "manifest_sha256": OWNER_MANIFEST,
        **{name: {"source": key(OWNER_RAW / name), **identity(WORK / name)} for name in ("run_owned.py", "owned_guard.py")}},
    "private_tmpdir": {"path": str(WORK / "temporary"), "created_only_at_execution": True,
        "mode": "0700", "owner_uid": os.getuid(), "preserve_natural_residuals": True}, "pins": pins,
})
print(json.dumps({"prepared_without_execution": True, "names": 9, "core": 8, "render": 1, "pins": len(pins),
                  "preparation": identity(WORK / "preparation.json"), "command": ["python3", "-B", key(WORK / "run_radiation.py")]}, indent=2))
