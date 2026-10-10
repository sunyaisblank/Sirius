"""Preserve the exact inputs of the one executed b728 wide software observer."""
from pathlib import Path
import hashlib, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent

def document(path):
    return json.loads(path.read_text())

def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

head = document(WORK / "source.json")["source_revision"]
assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == head
assert not subprocess.check_output(["git", "status", "--porcelain"])
assert document(WORK / "independent-source-stream-review.json")["pass"]
case = "software-wide-observer"
folder = WORK / case
owner = document(folder / "owner.json")
assert owner["passed"] and owner["child_exit"] == 0
assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
assert owner["stop_reason"] is None and not owner["cleanup_errors"]
assert not owner["remaining_owned_processes"]
assert all(v["absent"] for v in owner["observed_births"])
assert owner["source"]["revision"] == head and owner["source"]["status"] == ""
assert document(folder / "accepted.json")["pass"]
inputs = document(folder / "inputs-before.json")
assert inputs == document(folder / "inputs-after.json") and len(inputs) == 575
tree = {}
for row in subprocess.check_output(["git", "ls-tree", "-r", "-z", head]).split(b"\0"):
    if row:
        info, path = row.split(b"\t", 1)
        mode, kind, oid = info.decode().split()
        assert kind == "blob"
        tree[path.decode()] = oid

DEST = ROOT / "attestations/review" / head[:7] / "fma-dense-dp-software-observer"
assert not DEST.exists()
selected = {}

def select(path, target):
    assert path.is_file() and not path.is_symlink()
    assert target not in selected or identity(selected[target]) == identity(path)
    selected[target] = path

ledger = []
for relative, seal in sorted(inputs.items()):
    origin = ROOT / relative
    assert seal["resolved_path"] == str(origin.resolve())
    value = {k: seal[k] for k in ("bytes", "sha256")}
    if relative in tree:
        data = subprocess.check_output(["git", "cat-file", "blob", tree[relative]])
        assert {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()} == value
        recovery = {"kind": "immutable-Git-blob", "revision": head,
                    "blob": tree[relative], "path": relative}
    else:
        assert identity(origin) == value, relative
        target = "executed-inputs/" + value["sha256"]
        select(origin, target)
        recovery = {"kind": "local-archive-payload", "path": target}
    ledger.append({"origin": relative, **value, "recovery": recovery})
for path in folder.rglob("*"):
    assert not path.is_symlink()
    if path.is_file():
        select(path, "diagnostic/" + str(path.relative_to(WORK)))
for name in ("source.json", "actions.json", "baseline-arrays.json",
             "build-and-arrays.json", "linux-readonly-payloads.json",
             "independent-source-stream-review.json", "expected-readback-layout.json",
             "software-validator-rejected-metadata.json", "reader-finite-invalid-controls.json",
             "observer-result-comment.md", "observer-result-publication.json",
             "run_owned.py", "owned_guard.py", "bind_linux.py", "native_contracts.py",
             "wide-observation.patch", "protect_software.py"):
    select(WORK / name, "diagnostic/" + name)
payloads = []
for target, path in sorted(selected.items()):
    value = identity(path)
    output = DEST / target
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, output)
    assert identity(output) == value
    payloads.append({"path": target, "origin": str(path.relative_to(ROOT)), **value})
assert {str(p.relative_to(DEST)) for p in DEST.rglob("*") if p.is_file()} == {v["path"] for v in payloads}
for item in payloads:
    assert identity(DEST / item["path"]) == {k: item[k] for k in ("bytes", "sha256")}
for item in ledger:
    if item["recovery"]["kind"] == "local-archive-payload":
        assert identity(DEST / item["recovery"]["path"]) == {k: item[k] for k in ("bytes", "sha256")}
manifest = {
    "revision": head, "payloads": payloads, "payload_count": len(payloads),
    "payload_bytes": sum(v["bytes"] for v in payloads),
    "executed_cases": {case: {"source_revision": head, "passed": True, "inputs": ledger}},
    "scope": "Exact 575 input versions and original diagnostics of one executed wide software observer. "
             "Recovery from pushed source Git blobs and hashed local archive payloads. "
             "No recovery claim for overwritten d07 executables, sealed external software provider, "
             "native FMA, performance, frame or full qualification."
}
(DEST / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
receipt = {"archive": str(DEST.relative_to(ROOT)), "manifest": identity(DEST / "manifest.json"),
           "payload_count": len(payloads), "payload_bytes": manifest["payload_bytes"],
           "case_input_counts": {case: len(ledger)}, "fresh_byte_checks_pass": True}
(WORK / "software-protection.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
