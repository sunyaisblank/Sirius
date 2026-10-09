import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REVISION = os.environ.get("SIRIUS_DIAGNOSTIC_REVISION", "e0df788f67975047a62fa092a5a2943a75e26bba")
assert len(REVISION) == 40 and all(c in "0123456789abcdef" for c in REVISION)
IS_CURRENT = "SIRIUS_DIAGNOSTIC_REVISION" in os.environ
WORKFLOW_RUN = int(os.environ.get("SIRIUS_DIAGNOSTIC_WORKFLOW_RUN", "37103105221"))
EVIDENCE = ROOT / ("attestations/production-acceptance/" + REVISION[:7] if IS_CURRENT
                   else "attestations/retained-p1/e0df788")
BUNDLE = EVIDENCE / "windows-build"
RUNTIME = ROOT / "bin/windows-msvc" / ("production-acceptance-" + REVISION[:7] if IS_CURRENT
                                      else "retained-p1-e0df788")
POWERSHELL = "/mnt/c/WINDOWS/System32/WindowsPowerShell/v1.0/powershell.exe"


def windows_path(path):
    return subprocess.check_output(["wslpath", "-w", str(path)], text=True).strip()


def powershell(script, output):
    script_path = ROOT / "out/retained-qualification/native-check.ps1"
    script_path.write_text(
        "[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)\n"
        + script,
        encoding="utf-8",
    )
    with output.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            [POWERSHELL, "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
             "-File", windows_path(script_path)],
            stdout=stream, stderr=subprocess.STDOUT, cwd=ROOT,
        )
    print(f"native command exit={completed.returncode}; output={output.relative_to(ROOT)}", flush=True)
    return completed.returncode


def quote(value):
    return "'" + value.replace("'", "''") + "'"


def assemble():
    spec = importlib.util.spec_from_file_location("attestation", ROOT / "scripts/verify-attestation.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    document_path = BUNDLE / "windows-build.json"
    verifier.verify_path(document_path, ROOT)
    document = json.loads(document_path.read_text())
    assert document["source_revision"] == REVISION
    artifacts = document["artifacts"]
    copied = {}

    def copy(key, target):
        artifact = artifacts[key]
        source = BUNDLE / artifact["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        assert target.stat().st_size == artifact["bytes"] and digest == artifact["sha256"]
        copied[str(target.relative_to(RUNTIME))] = digest

    copy("qualification-sirius.bin", RUNTIME / "sirius.exe")
    copy("alignment_receipt.json", RUNTIME / "resources/model/alignment_receipt.json")
    for name, key in verifier.QUALIFICATION_PRODUCT_EVIDENCE.items():
        copy(key, RUNTIME / "resources" / verifier.QUALIFICATION_RUNTIME_RESOURCE_PATHS[name])
    for name in ("sirius_backend_tests", "sirius_render_tests"):
        copy(verifier.QUALIFICATION_TEST_EVIDENCE[name], RUNTIME / f"{name}.exe")
    for path in RUNTIME.glob("*.exe"):
        path.chmod(0o755)
    (EVIDENCE / "runtime-identity.json").write_text(json.dumps({
        "source_revision": REVISION,
        "workflow_run": WORKFLOW_RUN,
        "scope": "verified native compilation bundle; finite numerical diagnostics; no Mandatory/runtime/release promotion",
        "runtime_sha256": copied,
    }, indent=2) + "\n")
    inventory = EVIDENCE / "radeon-inventory.json"
    assert powershell(
        f"& {quote(windows_path(RUNTIME / 'sirius.exe'))} --json info system\nexit $LASTEXITCODE\n",
        inventory,
    ) == 0
    data = json.loads(inventory.read_text(encoding="utf-8-sig"))
    devices = data["backends"]["vulkan"]["devices"]
    physical = [d for d in devices if d["kind"] == "integrated" and "Radeon 780M" in d["name"]]
    assert len(physical) == 1
    print(json.dumps(physical[0], indent=2), flush=True)


def run(executable, test_filter, name, repeat="1", render_directory=None):
    repeat = int(repeat)
    assert repeat > 0
    data = json.loads((EVIDENCE / "radeon-inventory.json").read_text(encoding="utf-8-sig"))
    device = next(d for d in data["backends"]["vulkan"]["devices"] if "Radeon 780M" in d["name"])
    temporary = ROOT / "out/retained-qualification/native-temp"
    rendered = ROOT / "renders" / ("production-acceptance-" + REVISION[:7] if IS_CURRENT
                                     else "retained-p1-e0df788")
    if render_directory is not None:
        rendered = Path(render_directory)
        assert rendered.resolve().is_relative_to((ROOT / "renders").resolve())
    temporary.mkdir(parents=True, exist_ok=True)
    rendered.mkdir(parents=True, exist_ok=True)
    xml = EVIDENCE / f"{name}.xml"
    script = (
        f"$env:TEMP = {quote(windows_path(temporary))}\n"
        f"$env:TMP = $env:TEMP\n"
        f"$env:TMPDIR = {quote(windows_path(rendered))}\n"
        f"$env:SIRIUS_MEMORY_BUDGET_MB = '2048'\n"
        f"$env:SIRIUS_VULKAN_DEVICE = '{device['index']}'\n"
        f"& {quote(windows_path(RUNTIME / (executable + '.exe')))} "
        f"{quote('--gtest_filter=' + test_filter)} "
        f"{quote('--gtest_output=xml:' + windows_path(xml))} "
        f"{quote('--gtest_repeat=' + str(repeat))}\n"
        "exit $LASTEXITCODE\n"
    )
    sys.exit(powershell(script, EVIDENCE / f"{name}.log"))


if __name__ == "__main__":
    if sys.argv[1] == "assemble":
        assemble()
    else:
        run(*sys.argv[1:])
