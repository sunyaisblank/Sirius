#!/usr/bin/env python3
"""Verify complete CI build bundles against immutable Git source authority."""
import hashlib
import importlib.util
import json
import stat
import subprocess
import sys
import zipfile
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
WORK = Path(__file__).resolve().parent
REVISION = "5b5fad91c08c0533e29c4f5c7ae249725162e405"
RUN = 38577826600
PINNED = {
    "macos": (11901712901, 57718944, "fab6f99bd9a11b2e042b5d5f32d52e6cc25937e2aadcfe2077f35a4e4d9bb32a"),
    "windows": (11901703386, 56546068, "24b09db18d6989856e438c6c3dcfce941b5d87463dca7812ed2be199f7f093db"),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, timeout=30)


def save_json(path, data):
    with path.open("x") as stream:
        stream.write(json.dumps(data, indent=2, sort_keys=True) + "\n")


def main():
    require(not (WORK / "verification.json").exists(), "one-shot verification already attempted")
    before = git("rev-parse", "HEAD")
    require(not git("status", "--porcelain"), "runner checkout must be clean")
    require(git("rev-parse", REVISION + "^{commit}").decode().strip() == REVISION,
            "immutable source revision is absent")
    authority = WORK / "source-authority"
    authority.mkdir()
    source_records = {}
    for name in ("scripts/verify-attestation.py", "scripts/verify-build-gate.py",
                 "tests/operating_model.json", ".github/workflows/ci.yml"):
        payload = git("show", REVISION + ":" + name)
        target = authority / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        source_records[name] = {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}
        if name.startswith("scripts/"):
            require((ROOT / name).read_bytes() == payload, "runner verifier differs: " + name)
    model = (authority / "tests/operating_model.json").read_bytes()
    require(model.endswith(b"\n") and b"\r" not in model, "model is not canonical LF")
    decoded = json.loads(model)
    require(decoded["schema_version"] == 3 and decoded["method"] == "adversarial-claim-ledger",
            "unsupported operating model")
    model_sha = hashlib.sha256(model).hexdigest()
    path = authority / "scripts/verify-attestation.py"
    spec = importlib.util.spec_from_file_location("immutable_build_attestation_verifier", path)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    run_raw = subprocess.check_output([
        "gh", "api", f"repos/sunyaisblank/Sirius/actions/runs/{RUN}"], cwd=ROOT, timeout=60)
    (WORK / "run-api.json").write_bytes(run_raw)
    run = json.loads(run_raw)
    require(run["head_sha"] == REVISION and run["id"] == RUN, "CI run source mismatch")
    results = {}
    for role, (artifact_id, length, digest) in PINNED.items():
        raw = subprocess.check_output([
            "gh", "api", f"repos/sunyaisblank/Sirius/actions/artifacts/{artifact_id}"],
            cwd=ROOT, timeout=60)
        (WORK / (role + "-artifact-api.json")).write_bytes(raw)
        metadata = json.loads(raw)
        require(metadata["id"] == artifact_id and metadata["name"] == f"sirius-{role}-native-attestation"
                and metadata["size_in_bytes"] == length and not metadata["expired"]
                and metadata["workflow_run"]["id"] == RUN
                and metadata["workflow_run"]["head_sha"] == REVISION,
                "artifact API provenance mismatch: " + role)
        archive = WORK / (role + ".zip")
        payload = archive.read_bytes()
        require(len(payload) == length and hashlib.sha256(payload).hexdigest() == digest,
                "downloaded ZIP identity mismatch: " + role)
        require(metadata.get("digest") == "sha256:" + digest, "API digest mismatch: " + role)
        destination = WORK / role
        destination.mkdir()
        members = {}
        with zipfile.ZipFile(archive) as bundle:
            entries = bundle.infolist()
            require(len(entries) == 40 and len({i.filename for i in entries}) == 40,
                    "unexpected ZIP member inventory")
            require(sum(i.file_size for i in entries) < 256 * 1024 * 1024, "ZIP exceeds extraction bound")
            for info in entries:
                name = info.filename
                require(Path(name).name == name and name not in (".", "..")
                        and "/" not in name and "\\" not in name and ":" not in name
                        and stat.S_ISREG(info.external_attr >> 16)
                        and not info.is_dir() and not (info.flag_bits & 1),
                        "unsafe ZIP member: " + name)
                data = bundle.read(info)  # ZIP CRC is checked by the standard library.
                require(len(data) == info.file_size, "ZIP member size mismatch")
                with (destination / name).open("xb") as stream:
                    stream.write(data)
                members[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        location = destination / (role + "-build.json")
        document = json.loads(location.read_bytes())
        require(document["domains"] == [role + "-native-build"], "unexpected domain scope")
        require(set(members) == {location.name} | {r.get("path", n) for n, r in document["artifacts"].items()},
                "unclaimed ZIP members")
        verifier.verify_document_against_authority(document, location, REVISION, model_sha)
        results[role] = {"verified": True, "domains": document["domains"],
                         "claims": document["claims"], "members": members,
                         "archive": {"bytes": length, "sha256": digest}, "artifact_id": artifact_id}
        print(json.dumps({"role": role, "verified": True, "domains": document["domains"],
                          "test_report": document["claims"]["test_report"]}), flush=True)
    require(git("rev-parse", "HEAD") == before and not git("status", "--porcelain"),
            "runner checkout changed")
    save_json(WORK / "verification.json", {
        "passed": True, "source_revision": REVISION, "operating_model_sha256": model_sha,
        "verifier_api": "verify_document_against_authority",
        "source_authority": "exact immutable Git blobs; no runtime claim or worktree substitution",
        "source_files": source_records, "runner_head": before.decode().strip(),
        "ci_run_id": RUN, "results": results,
        "external_domains_admitted": 2, "external_domains_required": 8,
        "new_runtime_domain_admitted": False, "scientific_or_release_qualified": False,
    })


if __name__ == "__main__":
    main()
