"""Root-only single-action launcher; it replaces itself with the prior owner."""
from pathlib import Path
import os
import sys

sys.dont_write_bytecode = True
from check import preflight, validate

work = Path(__file__).resolve().parent
assert len(sys.argv) == 2
case = sys.argv[1]
assert case in {"inventory", "cie", "xyz"}
preflight(case)
assert not (work / case).exists(), "No retry or overwrite of a previous action"
assert not (work / ("accepted-" + case + ".json")).exists()
if case == "cie":
    validate("inventory")
if case == "xyz":
    validate("cie")
os.execv(sys.executable, [sys.executable, str(work / "run_owned.py"), case])
