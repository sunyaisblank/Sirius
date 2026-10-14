"""Preserve an unexecuted proposal and change only its outer time allowance."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import shutil
import sys

sys.dont_write_bytecode = True
WORK = Path(__file__).resolve().parent
ROOT = WORK.parents[1]
CASE = "original-doppler-raster-control"

def identity(path):
    raw = Path(path).read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}

def key(path):
    return str(path.relative_to(ROOT))

def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")

assert not any((WORK / name).exists() for name in (CASE, "result.json", "temporary"))
prepared = json.loads((WORK / "preparation.json").read_text())
actions = json.loads((WORK / "actions.json").read_text())
assert set(actions) == {CASE} and actions[CASE]["limits"] == {
    "seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}
assert prepared["executed"] is False
for name, pin in prepared["pins"].items():
    assert identity(ROOT / name) == pin, name

historical = WORK / "historical/proposal-300s"
assert not historical.exists()
historical.mkdir(parents=True)
records = {}
for name in ("prepare.py", "preparation.json", "actions.json", "run_doppler.py"):
    source, target = WORK / name, historical / name
    shutil.copy2(source, target)
    assert identity(source) == identity(target)
    records[name] = {"original_path": key(source), "preserved_path": key(target), **identity(target)}
write(historical / "preservation.json", {"unexecuted_proposal": True,
    "original_limits": actions[CASE]["limits"], "records": records,
    "relocation_scope": "The four changed proposal files are preserved byte-exact; every other original frozen input remains at its original unchanged location. No original Doppler attempt exists."})

wrapper = (WORK / "run_doppler.py").read_text()
old = 'LIMITS = {"seconds": 300, "rss_kib": 4194304, "reserve_kib": 2097152}'
new = old.replace('"seconds": 300', '"seconds": 900')
assert wrapper.count(old) == 1
(WORK / "run_doppler.py").write_text(wrapper.replace(old, new))

recipe = (WORK / "prepare.py").read_text()
assert recipe.count("'seconds': 300") == 1 and recipe.count("TIMEOUT;300s") == 1
recipe = recipe.replace("'seconds': 300", "'seconds': 900").replace("TIMEOUT;300s", "TIMEOUT;900s")
marker = "wrapper = wrapper_path.read_text()\n"
assert recipe.count(marker) == 1
recipe = recipe.replace(marker, marker +
    "assert wrapper.count('LIMITS = {\\\"seconds\\\": 300,') == 1\n" +
    "wrapper = wrapper.replace('LIMITS = {\\\"seconds\\\": 300,', 'LIMITS = {\\\"seconds\\\": 900,', 1)\n")
(WORK / "prepare.py").write_text(recipe)

prepared["reference"]["cost_scope"]["scope"] = (
    "Original registration has no TIMEOUT;900s is an operational admission ceiling only, "
    "not an original scientific bound or a completion forecast. Both complete120x68 passes "
    "and all original assertions remain unchanged; timeout means incomplete.")
prepared["outer_allowance_revision"] = {
    "previous_proposal_seconds": 300, "new_operational_ceiling_seconds": 900,
    "original_attempt_existed": False,
    "rationale": "Original source has16320 rays,20000 trace-attempt cap per ray and up to16 workers. Root reports the original577-trace/5000-step edge fixture took296.75s and576-search fixtures180–248s. These observations do not forecast Doppler completion;300s had no demonstrated adequacy.900s changes only operational admission, with no numerical/assertion/qualification/runtime-gate change or completion promise.",
    "historical_proposal": records,
    "preservation_receipt": {"path": key(historical / "preservation.json"), **identity(historical / "preservation.json")},
    "amendment_source": {"path": key(Path(__file__)), **identity(Path(__file__))}}
prepared["prepared_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
pins = prepared["pins"]
for path in [WORK / "prepare.py", WORK / "run_doppler.py", Path(__file__),
             *[historical / name for name in (*records, "preservation.json")]]:
    pins[key(path)] = identity(path)
actions[CASE]["limits"]["seconds"] = 900
actions[CASE]["frozen_inputs"] = list(dict.fromkeys([*pins, key(WORK / "preparation.json")]))
write(WORK / "actions.json", actions)
pins[key(WORK / "actions.json")] = identity(WORK / "actions.json")
write(WORK / "preparation.json", prepared)
for name in ("prepare.py", "run_doppler.py", "amend_outer_allowance.py"):
    ast.parse((WORK / name).read_text())
assert not any((WORK / name).exists() for name in (CASE, "result.json", "temporary"))
print(json.dumps({"prepared_without_test_execution": True, "pins": len(pins),
    "limits": actions[CASE]["limits"], "preparation": identity(WORK / "preparation.json"),
    "actions": identity(WORK / "actions.json"), "wrapper": identity(WORK / "run_doppler.py")}, indent=2))
