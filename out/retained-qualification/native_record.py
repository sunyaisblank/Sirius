"""Disposable identity wrapper for finite native Radeon diagnostics."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
revision = os.environ["SIRIUS_DIAGNOSTIC_REVISION"]
workflow = int(os.environ["SIRIUS_DIAGNOSTIC_WORKFLOW_RUN"])
executable, test_filter, name, scope = sys.argv[1:5]
render_directory = sys.argv[5] if len(sys.argv) > 5 else None
evidence = ROOT / "attestations/production-acceptance" / revision[:7]
runtime = ROOT / "bin/windows-msvc" / ("production-acceptance-" + revision[:7])
manifest = json.loads((evidence / "runtime-identity.json").read_text())
assert manifest["source_revision"] == revision and manifest["workflow_run"] == workflow

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify_runtime():
    for relative, expected in manifest["runtime_sha256"].items():
        assert digest(runtime / relative) == expected, relative

verify_runtime()
for suffix in (".log", ".xml", "-identity.json"):
    assert not (evidence / (name + suffix)).exists(), "refuse to overwrite previous evidence"
arguments = [sys.executable, str(ROOT / "out/retained-qualification/native_checks.py"),
             executable, test_filter, name, "1"]
if render_directory:
    arguments.append(render_directory)
started_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
started = time.monotonic()
completed = subprocess.run(arguments, cwd=ROOT)
elapsed = time.monotonic() - started
verify_runtime()
record = {
    "source_revision": revision, "workflow_run": workflow, "scope": scope,
    "binary_sha256": manifest["runtime_sha256"][executable + ".exe"],
    "filter": test_filter, "started_at_utc": started_at,
    "exit_code": completed.returncode, "harness_wall_seconds": elapsed,
}
for suffix, key in ((".log", "log_sha256"), (".xml", "xml_sha256")):
    path = evidence / (name + suffix)
    if path.exists():
        record[key] = digest(path)
xml = evidence / (name + ".xml")
if xml.exists():
    estate = ET.parse(xml).getroot()
    record["estate"] = estate.attrib
    record["skipped_testcases"] = len(estate.findall(".//skipped"))
if render_directory:
    rendered = Path(render_directory).resolve()
    assert rendered.is_relative_to((ROOT / "renders").resolve())
    record["render_directory"] = str(rendered.relative_to(ROOT))
    record["outputs"] = {
        str(path.relative_to(ROOT)): {"bytes": path.stat().st_size, "sha256": digest(path)}
        for path in sorted(rendered.glob("*.exr"))
    }
(evidence / (name + "-identity.json")).write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2), flush=True)
sys.exit(completed.returncode)
