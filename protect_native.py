"""Preserve native trial outcomes and every original execution input version."""
from pathlib import Path
import hashlib, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent

def read(path): return json.loads(path.read_text())
def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
def key(value): return value["bytes"], value["sha256"]

bindings = read(WORK / "execution-bindings.json")
head, baseline = bindings["candidate_revision"], bindings["baseline_revision"]
assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == head
assert not subprocess.check_output(["git", "status", "--porcelain"])
disposition = read(WORK / "native-disposition-reviewed.json")
assert disposition["source_revision"] == head and disposition["batch_terminal"]
DEST = ROOT / "attestations/native-vulkan" / head[:7] / "wide-dense-dp-trial"
assert not DEST.exists()
recoveries, payloads, trees = {}, {}, {}

software = ROOT / read(WORK / "software-protection.json")["archive"]
for record in read(software / "manifest.json")["payloads"]:
    path = software / record["path"]
    assert identity(path) == {k: record[k] for k in ("bytes", "sha256")}
    recoveries[key(record)] = {"kind": "retained-software-payload", "path": str(path.relative_to(ROOT))}
for revision in (baseline, head):
    for path in (ROOT / "attestations/native-build" / revision[:7]).rglob("*"):
        assert not path.is_symlink()
        if path.is_file():
            recoveries[key(identity(path))] = {"kind": "retained-native-build-payload", "path": str(path.relative_to(ROOT))}
    tree = {}
    for row in subprocess.check_output(["git", "ls-tree", "-r", "-z", revision]).split(b"\0"):
        if row:
            info, path = row.split(b"\t", 1)
            mode, kind, oid = info.decode().split()
            assert kind == "blob"
            tree[path.decode()] = oid
    trees[revision] = tree

prior = read(ROOT / "attestations/native-vulkan/95c36cb/original-fp64-observation/manifest.json")
wanted_provider = {key(v) for origin, v in bindings["files"].items() if origin.startswith("/mnt/c/")}
for record in prior["executed_input_ledger"]:
    if record["origin"].startswith("/mnt/c/") and key(record) in wanted_provider:
        recovery = record["recovery"]
        assert recovery["kind"] != "immutable-Git-blob"
        assert identity(ROOT / recovery["path"]) == {k: record[k] for k in ("bytes", "sha256")}
        recoveries[key(record)] = {"kind": "retained-provider-payload", "path": recovery["path"]}

tools = read(WORK / "tools-preservation.json")
assert tools["pushed"] and tools["remote_peeled_commit_verified"]
for record in tools["files"]:
    data = subprocess.check_output(["git", "cat-file", "blob", record["blob"]])
    assert {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()} == {k: record[k] for k in ("bytes", "sha256")}
    recoveries[key(record)] = {"kind": "immutable-Git-blob", "revision": tools["commit"], "blob": record["blob"], "path": record["path"]}

def preserve(path, value):
    assert path.is_file() and not path.is_symlink() and identity(path) == value, str(path)
    if key(value) in recoveries: return recoveries[key(value)]
    target = "payloads/" + value["sha256"]
    output = DEST / target
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, output)
    assert identity(output) == value
    payloads[target] = {"path": target, **value}
    recovery = {"kind": "local-archive-payload", "path": str(output.relative_to(ROOT))}
    recoveries[key(value)] = recovery
    return recovery

ledger = []
for origin, value in sorted(bindings["files"].items()):
    path = ROOT / origin
    assert identity(path) == value, origin
    revision, source_path = head, origin
    for candidate in (baseline, head):
        prefix = "bin/windows-msvc/native-" + candidate[:7] + "/recorded-source/"
        if origin.startswith(prefix): revision, source_path = candidate, origin[len(prefix):]
    if source_path in trees[revision]:
        blob = trees[revision][source_path]
        assert subprocess.check_output(["git", "hash-object", "--no-filters", str(path)], text=True).strip() == blob
        recovery = {"kind": "immutable-Git-blob", "revision": revision, "blob": blob, "path": source_path}
    else:
        recovery = preserve(path, value)
    ledger.append({"origin": origin, **value, "recovery": recovery})
diagnostics = []
for path in sorted(WORK.rglob("*")):
    assert not path.is_symlink()
    if path.is_file():
        value = identity(path)
        diagnostics.append({"origin": str(path.relative_to(ROOT)), **value, "recovery": preserve(path, value)})
for record in (*ledger, *diagnostics):
    recovery = record["recovery"]
    expected = {k: record[k] for k in ("bytes", "sha256")}
    if recovery["kind"] == "immutable-Git-blob":
        data = subprocess.check_output(["git", "cat-file", "blob", recovery["blob"]])
        assert {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()} == expected
    else:
        assert identity(ROOT / recovery["path"]) == expected
    assert identity(ROOT / record["origin"]) == expected
assert {str(p.relative_to(DEST)) for p in DEST.rglob("*") if p.is_file()} == set(payloads)
manifest = {
    "source_revision": head, "baseline_revision": baseline, "disposition": disposition,
    "executed_input_ledger": ledger, "diagnostics": diagnostics,
    "payloads": list(payloads.values()), "payload_count": len(payloads),
    "payload_bytes": sum(v["bytes"] for v in payloads.values()),
    "scope": "Exact frozen trial inputs and original raw outcomes, including rejected outcomes if present. "
             "Separate Windows/WSL clocks remain uncalibrated; ownership visibility is bounded/accessible. "
             "No claim of whole scientific estate, frame, cold state or full release qualification; "
             "no recovery claim for overwritten d07 executables and no cleanup clearance."
}
(DEST / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
receipt = {"archive": str(DEST.relative_to(ROOT)), "manifest": identity(DEST / "manifest.json"),
           "payload_count": len(payloads), "payload_bytes": manifest["payload_bytes"],
           "executed_input_count": len(ledger), "diagnostic_count": len(diagnostics),
           "all_origin_and_recovery_bytes_verified": True, "cleanup_authorized": False}
(WORK / "native-protection.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
