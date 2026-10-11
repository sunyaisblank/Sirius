"""Preserve exact inputs of the three restored-source finite checks."""
from pathlib import Path
import hashlib, json, shutil, subprocess, xml.etree.ElementTree as ET

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
CHECKS = WORK / "restoration"

def document(path):
    return json.loads(path.read_text())

def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

head = document(CHECKS / "source.json")["source_revision"]
assert subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip() == head
assert not subprocess.check_output(["git", "status", "--porcelain"])
assert document(WORK / "restoration-source-independent-review.json")["pass"]
build = document(CHECKS / "build-and-arrays.json")
assert build["source_revision"] == head and build["original40_complete_arrays_and_whole_header_byte_exact95"]
assert len(build["arrays"]) == 40 and build["declared_spv_count"] == 33
bound = document(CHECKS / "linux-readonly-payloads.json")
assert bound["source_revision"] == head and bound["pass_"]
assert len(bound["actual_consumers"]) == 4
assert all(len(v["arrays"]) == 40 for v in bound["actual_consumers"].values())
xml = ET.parse(CHECKS / "admission/gtest.xml").getroot()
assert xml.get("tests") == "2" and all(int(xml.get(k, "0")) == 0 for k in ("failures", "errors", "disabled"))
expected = {"RetainedComputeAdmission.ArithmeticRefusalPrecedesKernelLoading",
            "RetainedComputeAdmission.FmaSelectsOnlyNativeWideProductsAndPreservesAllocation"}
cases_xml = xml.findall(".//testcase")
assert len(cases_xml) == 2 and {v.get("classname") + "." + v.get("name") for v in cases_xml} == expected
assert all(v.get("status") == "run" and v.get("result") == "completed" and v.find("skipped") is None for v in cases_xml)
cases = ("configure", "build", "admission")
ledgers = {}
tree = {}
for row in subprocess.check_output(["git", "ls-tree", "-r", "-z", head]).split(b"\0"):
    if row:
        info, path = row.split(b"\t", 1)
        mode, kind, oid = info.decode().split()
        assert kind == "blob"
        tree[path.decode()] = oid

DEST = ROOT / "attestations/review" / head[:7] / "rejected-upmultiply-restoration-source-checks"
assert not DEST.exists()
selected = {}

def select(path, target):
    assert path.is_file() and not path.is_symlink()
    assert target not in selected or identity(selected[target]) == identity(path)
    selected[target] = path

for case in cases:
    folder = CHECKS / case
    owner = document(folder / "owner.json")
    assert owner["passed"] and owner["child_exit"] == 0
    assert owner["source_unchanged"] and owner["whole_inputs_unchanged"]
    assert owner["stop_reason"] is None and not owner["cleanup_errors"]
    assert not owner["remaining_owned_processes"]
    assert all(v["absent"] for v in owner["observed_births"])
    assert owner["source"] == {"revision": head, "status": ""}
    inputs = document(folder / "inputs-before.json")
    assert inputs == document(folder / "inputs-after.json")
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
    ledgers[case] = {"source_revision": head, "passed": True, "inputs": ledger}
for path in CHECKS.rglob("*"):
    assert not path.is_symlink()
    if path.is_file():
        select(path, "diagnostic/restoration/" + str(path.relative_to(CHECKS)))
for name in ("restoration-source-contract.json", "restoration-source-independent-review.json",
             "protect_restoration.py", "bind_restoration.py"):
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
for case in cases:
    for item in ledgers[case]["inputs"]:
        if item["recovery"]["kind"] == "local-archive-payload":
            assert identity(DEST / item["recovery"]["path"]) == {k: item[k] for k in ("bytes", "sha256")}
manifest = {
    "revision": head, "payloads": payloads, "payload_count": len(payloads),
    "payload_bytes": sum(v["bytes"] for v in payloads),
    "executed_cases": ledgers,
    "scope": "Recoverable restored-source strict configure/build and two original host-model cases. "
             "All40 original payloads/header byte exact95,160 actual readonlyELFjoins. "
             "No current Vulkan/native/scientific/frame or full qualification claim."

}
(DEST / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
receipt = {"archive": str(DEST.relative_to(ROOT)), "manifest": identity(DEST / "manifest.json"),
           "payload_count": len(payloads), "payload_bytes": manifest["payload_bytes"],
           "case_input_counts": {case: len(ledgers[case]["inputs"]) for case in cases}, "fresh_byte_checks_pass": True}
(WORK / "restoration-protection.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt))
